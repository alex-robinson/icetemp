module icetemp_mali 
    ! Wrapping the MALI icetemp solution for yelmo 

    use defs, only : prec, pi, g, sec_year, T0, rho_ice, rho_sw, rho_w 
    use solver_tridiagonal, only : tridiag 
    use thermodynamics, only : calc_advec_horizontal_column_sico1, calc_advec_horizontal_column

    implicit none
    
    ! Impose parameter choice here for testing conductive bedrock or not
    logical,    parameter :: conductive_bedrock = .FALSE. 
    
    private
    !public :: calc_icetemp_mali_3D_up
    public :: calc_icetemp_mali_column_up
     




contains 

    subroutine calc_icetemp_mali_column_up(T_ice,bmb,is_float, H_ice, T_srf, advecxy, ux, uy, uz, dzsdx, dzsdy, dzsrfdt, &
                                      dHicedx, dHicedy, dHicedt, cp, kt, Q_strn, Q_b, T_pmp, mb_net, Q_geo, sigma, dx, dt)
        ! Calculate input variables to column temperature solver,
        ! and reverse index order (1:base,nz:surface) => (1:surface,nz:base)

        implicit none 

        real(prec), intent(INOUT) :: T_ice(:)                           ! The ice temperature [K]
        real(prec), intent(OUT)   :: bmb                               ! The basal mass balance (negative melt) at the bottom at time time + dt if T_pmp is reached [J m^-2 y^-1]
        logical,    intent(IN)  :: is_float                            ! Grounded fraction 
        real(prec), intent(IN)  :: H_ice                             ! The ice thickness [m]
        real(prec), intent(IN)  :: T_srf                             ! The ice surface temperature at time step time + dt [K]
        real(prec), intent(IN)  :: advecxy(:)                         ! The 3D velocity field in the x-direction [m y^-1]
        real(prec), intent(IN)  :: ux(:)                              ! The 3D velocity field in the x-direction [m y^-1]
        real(prec), intent(IN)  :: uy(:)                              ! The 3D velocity field in the y-direction [m y^-1]
        real(prec), intent(IN)  :: uz(:)                              ! The 3D velocity field in the zeta-direction [m y^-1]
        real(prec), intent(IN)  :: dzsdx                             ! The surface gradient in the x-direction [m m^-1]
        real(prec), intent(IN)  :: dzsdy                             ! The surface gradient in the y-direction [m m^-1]
        real(prec), intent(IN)  :: dzsrfdt                           ! The surface gradient in time [m a-1]
        real(prec), intent(IN)  :: dHicedx                           ! The ice thickness gradient in the x-direction [m m^-1]
        real(prec), intent(IN)  :: dHicedy                           ! The ice thickness gradient in the y-direction [m m^-1]
        real(prec), intent(IN)  :: dHicedt                           ! The ice thickness gradient in time [m a-1]
        real(prec), intent(IN)  :: cp(:)                              ! The specific heat capacity of ice at each x,y,zeta point [J kg^-1 K^-1]
        real(prec), intent(IN)  :: kt(:)                              ! The conductivity of ice at each x,y,zeta point [J m^-1 K^-1 y^-1]
        real(prec), intent(IN)  :: Q_strn(:)                          ! Internal strain heating [K a-1]
        real(prec), intent(IN)  :: Q_b                               ! The heat flux at the ice bottom [J m^-2 y^-1]
        real(prec), intent(IN)  :: T_pmp(:)                           ! The pressure melting point temperature for each depth and for all grid points [K]
        real(prec), intent(IN)  :: mb_net                            ! Net basal and surface mass balance [meter ice equivalent per year]
        real(prec), intent(IN)  :: Q_geo                             ! Geothermal heat flux [1e-3 J m^-2 s^-1]
        real(prec), intent(IN)  :: sigma(:)                               ! Vertical height axis (0:1) 
        real(prec), intent(IN)  :: dx 
        real(prec), intent(IN)  :: dt 

        ! Column variables for interfacing with heiko's routine 
        real(prec), allocatable  :: rev_Ti(:)                             ! The ice temperature at the previous time step time [K]
        real(prec), allocatable  :: rev_U(:)                              ! The 3D velocity field in the x-direction [m y^-1]
        real(prec), allocatable  :: rev_V(:)                              ! The 3D velocity field in the y-direction [m y^-1]
        real(prec), allocatable  :: rev_W(:)                              ! The 3D velocity field in the zeta-direction [m y^-1]
        real(prec), allocatable  :: rev_Q_strn(:)                         ! Reversed strain heating field [K y^-1]
        real(prec), allocatable  :: rev_dV_dzeta(:)                       ! The zeta-derivative of the 3D velocity field in the y-direction [y^-1]
        real(prec), allocatable  :: rev_Cpi(:)                            ! The specific heat capacity of ice at each x,y,zeta point [J kg^-1 K^-1]
        real(prec), allocatable  :: rev_Ki(:)                             ! The conductivity of ice at each x,y,zeta point [J m^-1 K^-1 y^-1]
        real(prec), allocatable  :: rev_Ti_pmp(:)                         ! The pressure melting point temperature for each depth and for all grid points [K]
        real(prec), allocatable  :: rev_advecxy(:)                        ! The xy horizontal advection of temperature contribution
        real(prec), allocatable  :: zeta(:)                               ! The xy horizontal advection of temperature contribution
        
        ! Output variables:
        real(prec), allocatable  :: rev_Ti_new(:)                          ! The new ice temperature at time time + dt [K]. The surface temperature T(1,:,:) is calculated with ice_surface_temperature()
        
        ! Local variables
        integer    :: i, j, k, nx, ny, nz  
        real(prec) :: bottom_melt 
        real(prec) :: Q_geo_now 
        real(prec) :: ghf_conv, ghf_now 

        real(prec), allocatable :: zeta_t(:)
        real(prec), allocatable :: zeta_x(:)
        real(prec), allocatable :: zeta_y(:)
        real(prec) :: zeta_z

        !type(zeta_helper_type)  :: Cz                    ! Zeta helpers 
        
        ghf_conv = 1e-3*sec_year   ! [mW m-2] => [J m-2 a-1]

        nz = size(T_ice,1)

        allocate(zeta(nz)) 

        allocate(zeta_t(nz))
        allocate(zeta_x(nz))
        allocate(zeta_y(nz))

        allocate(rev_Ti(nz))
        allocate(rev_U(nz))
        allocate(rev_V(nz))
        allocate(rev_W(nz))
        allocate(rev_Q_strn(nz))
        allocate(rev_Cpi(nz))
        allocate(rev_Ki(nz))
        allocate(rev_Ti_pmp(nz))
        allocate(rev_advecxy(nz))
        allocate(rev_Ti_new(nz))

        ! Calculate depth axis from height axis (0:1) => (1:0) 
        ! Store reversed column variables
        do k = 1, nz 
            zeta(nz+1-k) = 1.0 - sigma(k)
        end do  

!         ! Calculate the zeta derivatives needed for thermodynamic solver
!         call calculate_zeta_derivatives(zeta_t,zeta_x,zeta_y,zeta_z, H_ice, &
!                             dHicedt, dHicedx, dHicedy, dzsrfdt, dzsdx, dzsdy, zeta)

!         ! Calculate zeta helpers 
!         call calc_zeta_helpers(Cz,zeta,dx,dx)

        ! Get geothermal heat flux in proper units 
        Q_geo_now = Q_geo*ghf_conv 

        ! Store reversed column variables
        do k = 1, nz  
            rev_Ti(nz+1-k)       = T_ice(k)
            rev_U(nz+1-k)        = ux(k)
            rev_V(nz+1-k)        = uy(k) 
            rev_W(nz+1-k)        = uz(k)
            rev_Q_strn(nz+1-k)   = Q_strn(k)
            rev_Cpi(nz+1-k)      = cp(k)
            rev_Ki(nz+1-k)       = kt(k)
            rev_Ti_pmp(nz+1-k)   = T_pmp(k)
            rev_advecxy(nz+1-k)  = advecxy(k)
        end do  
     
!         call temperature_mali(zeta, Cz, zeta_t(:), zeta_x(:), zeta_y(:), zeta_z, &
!                         is_float, H_ice, T_srf, rev_Ti, rev_U, rev_V, rev_W, &
!                         rev_Q_strn, dzsdx, dzsdy, rev_Cpi, rev_Ki, &
!                         Q_geo_now, rev_Ti_pmp, rev_advecxy,mb_net, &
!                         rev_Ti_new,bottom_melt,dt) 

        ! Return output variables back to original vertical ordering
        bmb = -bottom_melt 
        do k = 1, nz 
            T_ice(k) = rev_Ti_new(nz+1-k)
        end do 

        return 

    end subroutine calc_icetemp_mali_column_up






end module icetemp_mali


