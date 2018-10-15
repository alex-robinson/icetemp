module icetemp_imau
! === from IMAU-ICE (Heiko Goelzer) === 
    
    use defs, only : prec, pi, g, rho_ice, sec_year, T0, rho_w

    use solver_tridiagonal, only : solve_tridiag
    
    use thermodynamics 

    implicit none

    type zeta_helper_type

        real(prec), allocatable :: a_zeta(:) 
        real(prec), allocatable :: a_zeta_minus(:) 
        real(prec), allocatable :: a_zetazeta(:) 
        real(prec), allocatable :: b_zeta(:) 
        real(prec), allocatable :: b_zeta_minus(:) 
        real(prec), allocatable :: b_zetazeta(:) 
        real(prec), allocatable :: c_zeta(:) 
        real(prec), allocatable :: c_zetazeta(:) 
        real(prec), allocatable :: z_zeta_minus(:) 

    end type 

    private 
    public :: calc_icetemp_imau_3D_up
    public :: calc_icetemp_imau_column_up
    
contains 
    
    subroutine calc_icetemp_imau_3D_up(T_ice,bmb,f_grnd, H_ice, T_srf, ux, uy, uz, dzsdx, dzsdy, dzsrfdt, &
                                      dHicedx, dHicedy, dHicedt, cp, kt, Q_strn, Q_b, T_pmp, smb, Q_geo, sigma, dx, dt)
        ! Calculate input variables to column temperature solver,
        ! and reverse index order (1:base,nz:surface) => (1:surface,nz:base)

        implicit none 

        real(prec), intent(OUT) :: T_ice(:,:,:)                           ! The ice temperature [K]
        real(prec), intent(OUT) :: bmb(:,:)                               ! The basal mass balance (negative melt) at the bottom at time time + dt if T_pmp is reached [J m^-2 y^-1]
        real(prec), intent(IN)  :: f_grnd(:,:)                            ! Grounded fraction 
        real(prec), intent(IN)  :: H_ice(:,:)                             ! The ice thickness [m]
        real(prec), intent(IN)  :: T_srf(:,:)                             ! The ice surface temperature at time step time + dt [K]
        real(prec), intent(IN)  :: ux(:,:,:)                              ! The 3D velocity field in the x-direction [m y^-1]
        real(prec), intent(IN)  :: uy(:,:,:)                              ! The 3D velocity field in the y-direction [m y^-1]
        real(prec), intent(IN)  :: uz(:,:,:)                              ! The 3D velocity field in the zeta-direction [m y^-1]
        real(prec), intent(IN)  :: dzsdx(:,:)                             ! The surface gradient in the x-direction [m m^-1]
        real(prec), intent(IN)  :: dzsdy(:,:)                             ! The surface gradient in the y-direction [m m^-1]
        real(prec), intent(IN)  :: dzsrfdt(:,:)                           ! The surface gradient in time [m a-1]
        real(prec), intent(IN)  :: dHicedx(:,:)                           ! The ice thickness gradient in the x-direction [m m^-1]
        real(prec), intent(IN)  :: dHicedy(:,:)                           ! The ice thickness gradient in the y-direction [m m^-1]
        real(prec), intent(IN)  :: dHicedt(:,:)                           ! The ice thickness gradient in time [m a-1]
        real(prec), intent(IN)  :: cp(:,:,:)                              ! The specific heat capacity of ice at each x,y,zeta point [J kg^-1 K^-1]
        real(prec), intent(IN)  :: kt(:,:,:)                              ! The conductivity of ice at each x,y,zeta point [J m^-1 K^-1 y^-1]
        real(prec), intent(IN)  :: Q_strn(:,:,:)                          ! Internal strain heating [K a-1]
        real(prec), intent(IN)  :: Q_b(:,:)                               ! The heat flux at the ice bottom [J m^-2 y^-1]
        real(prec), intent(IN)  :: T_pmp(:,:,:)                           ! The pressure melting point temperature for each depth and for all grid points [K]
        real(prec), intent(IN)  :: smb(:,:)                               ! Surface mass balance [meter ice equivalent per year]
        real(prec), intent(IN)  :: Q_geo(:,:)                             ! Geothermal heat flux [1e-3 J m^-2 s^-1]
        real(prec), intent(IN)  :: sigma(:)                               ! Vertical height axis (0:1) 
        real(prec), intent(IN)  :: dx 
        real(prec), intent(IN)  :: dt 

        ! Local variables
        integer    :: i, j, nx, ny, nz 
        logical    :: is_float 
        real(prec), allocatable :: advecxy(:) 

        nx = size(T_ice,1)
        ny = size(T_ice,2)
        nz = size(T_ice,3)

        allocate(advecxy(nz))

        do j = 3, ny-2 
        do i = 3, nx-2 
        
            ! Determine whether current point is floating or grounded  
            is_float = f_grnd(i,j) .eq. 0.0 

            ! Calculate horizontal advection 
            call calc_advec_horizontal_column(advecxy,T_ice,ux,uy,dx,i,j)

            call calc_icetemp_imau_column_up(T_ice(i,j,:),bmb(i,j),is_float,H_ice(i,j),T_srf(i,j),advecxy, &
                                             ux(i,j,:),uy(i,j,:),uz(i,j,:),dzsdx(i,j),dzsdy(i,j),dzsrfdt(i,j), &
                                             dHicedx(i,j),dHicedy(i,j),dHicedt(i,j),cp(i,j,:),kt(i,j,:), &
                                             Q_strn(i,j,:),Q_b(i,j),T_pmp(i,j,:),smb(i,j)+bmb(i,j),Q_geo(i,j),sigma,dx,dt)


        end do 
        end do 

        return 

    end subroutine calc_icetemp_imau_3D_up

    subroutine calc_icetemp_imau_column_up(T_ice,bmb,is_float, H_ice, T_srf, advecxy, ux, uy, uz, dzsdx, dzsdy, dzsrfdt, &
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

        type(zeta_helper_type)  :: Cz                    ! Zeta helpers 
        
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

        ! Calculate the zeta derivatives needed for thermodynamic solver
        call calculate_zeta_derivatives(zeta_t,zeta_x,zeta_y,zeta_z, H_ice, &
                            dHicedt, dHicedx, dHicedy, dzsrfdt, dzsdx, dzsdy, zeta)

        ! Calculate zeta helpers 
        call calc_zeta_helpers(Cz,zeta,dx,dx)

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
     
        call temperature_imau(zeta, Cz, zeta_t(:), zeta_x(:), zeta_y(:), zeta_z, &
                        is_float, H_ice, T_srf, rev_Ti, rev_U, rev_V, rev_W, &
                        rev_Q_strn, dzsdx, dzsdy, rev_Cpi, rev_Ki, &
                        Q_geo_now, rev_Ti_pmp, rev_advecxy,mb_net, &
                        rev_Ti_new,bottom_melt,dt) 

        ! Return output variables back to original vertical ordering
        bmb = -bottom_melt 
        do k = 1, nz 
            T_ice(k) = rev_Ti_new(nz+1-k)
        end do 

        return 

    end subroutine calc_icetemp_imau_column_up


    subroutine temperature_imau(zeta, Cz, zeta_t, zeta_x, zeta_y, zeta_z, &
                            is_float, Hi, Ts, Ti, U, V, W, Q_strn, dHs_dx, dHs_dy, Cpi, Ki, &
                            q_bottom, Ti_pmp, advecxy, mb_net, Ti_new, bottom_melt,dt)
        ! This subroutine solves the temperature field with the thermodynamical equation (X.11).

        ! Note: column indices are 1:surface, nz:base and 
        ! zeta is the sigma coordinate for depth, so zeta=0:surface and zeta=1:base 

        implicit none

        ! Input variables:
        real(prec), intent(IN)  :: zeta(:)                           ! Depth axis (0:surface, 1:depth)
        type(zeta_helper_type), intent(IN)  :: Cz                    ! Zeta helpers 
        real(prec), intent(IN)  :: zeta_t(:)                         ! zeta derivative t
        real(prec), intent(IN)  :: zeta_x(:)                         ! zeta derivative x
        real(prec), intent(IN)  :: zeta_y(:)                         ! zeta derivative y
        real(prec), intent(IN)  :: zeta_z                            ! zeta derivative z
        logical,    intent(IN)  :: is_float                          ! Mask floating (true) or grounded (false) ice (ie, is_float)
        real(prec), intent(IN)  :: Hi                                ! The ice thickness [m]
        real(prec), intent(IN)  :: Ts                                ! The ice surface temperature at time step time + dt [K]
        real(prec), intent(IN)  :: Ti(:)                             ! The ice temperature at the previous time step time [K]
        real(prec), intent(IN)  :: U(:)                              ! The 3D velocity field in the x-direction [m y^-1]
        real(prec), intent(IN)  :: V(:)                              ! The 3D velocity field in the y-direction [m y^-1]
        real(prec), intent(IN)  :: W(:)                              ! The 3D velocity field in the zeta-direction [m y^-1]
        real(prec), intent(IN)  :: Q_strn(:)                         ! The internal strain heating rate [K y^-1]
        real(prec), intent(IN)  :: dHs_dx                            ! The surface gradient in the x-direction [m m^-1]
        real(prec), intent(IN)  :: dHs_dy                            ! The surface gradient in the y-direction [m m^-1]
        real(prec), intent(IN)  :: Cpi(:)                            ! The specific heat capacity of ice at each x,y,zeta point [J kg^-1 K^-1]
        real(prec), intent(IN)  :: Ki(:)                             ! The conductivity of ice at each x,y,zeta point [J m^-1 K^-1 y^-1]
        real(prec), intent(IN)  :: q_bottom                          ! The heat flux at the ice bottom [J m^-2 y^-1]
        real(prec), intent(IN)  :: Ti_pmp(:)                         ! The pressure melting point temperature for each depth and for all grid points [K]
        real(prec), intent(IN)  :: advecxy(:)                        ! The pressure melting point temperature for each depth and for all grid points [K]
        real(prec), intent(IN)  :: mb_net                            ! Net surface + basal Mass Balance [meter ice equivalent per year]
        real(prec), intent(IN)  :: dt                                ! [a] Time step 

        ! Output variables:
        real(prec), intent(OUT) :: Ti_new(:)                         ! The new ice temperature at time time + dt [K]. The surface temperature T(1,:,:) is calculated with ice_surface_temperature()
        real(prec), intent(OUT) :: bottom_melt                       ! The melt at the bottom at time time + dt if Ti_pmp is reached [J m^-2 y^-1]

        ! Local variables:
        integer                 :: i, j, k, nz 
        real(prec)              :: internal_heating
        real(prec)              :: f1, f2, f3
        real(prec), allocatable :: alpha(:)
        real(prec), allocatable :: beta(:)
        real(prec), allocatable :: gamma(:)
        real(prec), allocatable :: delta(:)
        real(prec)              :: Cpi_dependent_part
        real(prec)              :: Delta_T
        integer                 :: additional_shelf_index
        real(prec), parameter :: maximum_bottom_melt = 0.5          ! [m/a]

        real(prec), parameter :: triple_point_of_water = 273.160
        real(prec), parameter :: seawater_temperature  = 271.150
        real(prec), parameter :: Hi_min                = 0.01       ! [m]

        logical    :: linear_temperature_at_shelf
        logical    :: fixed_Ti_at_shelf_bottom
        integer    :: choice_thermodynamics_scheme
        integer    :: number_of_potential_bottom_melt_layers
        real(prec) :: c_0_specific_heat
        real(prec) :: c_Delta_specific_heat
        real(prec) :: gravity_constant
        real(prec) :: registration_limit_low_Ti
        real(prec) :: latent_heat
        real(prec) :: ice_density 

        real(prec) :: G_b_t 

        linear_temperature_at_shelf             = .TRUE. 
        fixed_Ti_at_shelf_bottom                = .TRUE. 
        choice_thermodynamics_scheme            = 1
        number_of_potential_bottom_melt_layers  = 14 
        c_0_specific_heat                       = 2127.5  
        c_Delta_specific_heat                   = 7.253
        gravity_constant                        = 9.81
        registration_limit_low_Ti               = 200.0 
        latent_heat                             = 333500.0 
        ice_density                             = rho_ice 
        
        G_b_t = 1/dt 

        nz = size(Ti_new)

        allocate(alpha(nz))
        allocate(beta(nz))
        allocate(gamma(nz))
        allocate(delta(nz))
        
        alpha = 0.0 
        beta  = 0.0 
        gamma = 0.0 
        delta = 0.0 

        IF(Hi .le. Hi_min) THEN
            ! At grid points with minimal ice thickness the temperature field at all vertical layers will be initialized at Ts:
            Ti_new(:) = Ts
        ELSE IF(linear_temperature_at_shelf .AND. is_float) THEN
            ! This option initializes the temperature at shelf points linear between Ts and seawater temperature:
            Ti_new(:) = Ts + zeta(:) * (seawater_temperature - Ts)
        ELSE
            ! Ice surface boundary condition, see equations (X.31 - X.33):
            beta (1) = 1.0
            gamma(1) = 0.0
            delta(1) = Ts

            ! Loop over the vertical but without the surface (k=1) and the bottom (k=NZ):
            DO k = 2, nz-1
                SELECT CASE(choice_thermodynamics_scheme)
                    CASE(1)
                        Cpi_dependent_part = Cpi(k)
                    CASE(2)
                        Cpi_dependent_part = (2.0 * Cpi(k) - c_0_specific_heat &
                                        + c_Delta_specific_heat * triple_point_of_water)
                    CASE DEFAULT
                        STOP 'STOP: Invalid value for:  choice_thermodynamics_scheme_config'
                END SELECT

                IF(.not. is_float) THEN
                    internal_heating = Q_strn(k)
!                     internal_heating = ((- gravity_constant * zeta(k)) / Cpi_dependent_part) &
!                                 * ( dU_dzeta(k) * dHs_dx + dV_dzeta(k) * dHs_dy )
                ELSE
                    internal_heating = 0.0
                END IF

                ! See equations (X.20 - X.22):
                f1 = (Ki(k) / (ice_density * Cpi_dependent_part)) * zeta_z**2

                f2 = zeta_t(k) + zeta_x(k) * U(k) + zeta_y(k) * V(k) + zeta_z * W(k)

                f3 = internal_heating + advecxy(k) - G_b_t * Ti(k)

                ! See equations (X.26 - X.29):
                alpha(k) = f1 * Cz%a_zetazeta(k) - f2 * Cz%a_zeta(k)
                beta (k) = f1 * Cz%b_zetazeta(k) - f2 * Cz%b_zeta(k) - G_b_t
                gamma(k) = f1 * Cz%c_zetazeta(k) - f2 * Cz%c_zeta(k)
                delta(k) = f3
            END DO ! End k loop

            IF(fixed_Ti_at_shelf_bottom .AND. is_float) THEN
                ! With this option the seawater temperature is used as the bottom boundary condition (see equations X.61 - X.63):
                alpha(nz) = 0.0
                beta (nz) = 1.0
                delta(nz) = seawater_temperature
            ELSE
                ! The default boundary condition (second order gradient) at the ice bottom (see equations (X.42 - X.44)):
                alpha(nz) = Cz%z_zeta_minus(nz) * beta (nz-1) - Cz%a_zeta_minus(nz) * alpha(nz-1)
                beta (nz) = Cz%z_zeta_minus(nz) * gamma(nz-1) - Cz%b_zeta_minus(nz) * alpha(nz-1) 
                delta(nz) = Cz%z_zeta_minus(nz) * delta(nz-1) + (q_bottom / (zeta_z * Ki(nz))) * alpha(nz-1)
            END IF

            ! Solving the tridiagonal set of equations (X.25) which contains the implicit scheme in the vertical per horizontal grid point:
!             Ti_new(:) = tridiagonal_solve(alpha, beta, gamma, delta, 'thermodynamics_module [temperature]')
            call solve_tridiag(alpha,beta,gamma,delta,Ti_new,nz)

        END IF ! End IF Hi_min

        bottom_melt = 0.0
        
        IF(is_float) THEN
            additional_shelf_index = -1
        ELSE
            additional_shelf_index = 0
        END IF

        IF(Hi .le. Hi_min) THEN
            ! Do nothing, saving computational time in checking
        ELSE IF(linear_temperature_at_shelf .AND. is_float) THEN
            ! Do nothing, saving computational time in checking -- this ELSE IF-statement can be commented if the Ti_pmp check is desired instead.
        ELSE
            
            DO k = 2, nz
                !Delta_T =           Ti_new(k) - Ti_pmp(k)                                       ! original method: taking the lower layer value
                Delta_T = 0.50 * (Ti_new(k) - Ti_pmp(k) + Ti_new(k-1) - Ti_pmp(k-1))  ! new      method: taking the average value of the zeta layer

                ! If the new ice temperature exceeds the local pressure melting point (see equation (X.45)) due to the energy flux, it will be limited
                ! to this local pressure melting point and the energy excess is used for ice melt at the bottom:
                IF(Delta_T > 0.0) THEN

                    IF(nz - k + 1 <= number_of_potential_bottom_melt_layers) THEN
                        ! In all layers where the temperature exceeds the pressure melting point, the energy leavings are calculated from the temperature differences,
                        ! and the summation of them over the vertical layers is applied as bottom melt per (thermo) time step (see equation (X.48)):
                        !bottom_melt = bottom_melt + Delta_T * (zeta(k) - zeta(k-1)) * Hi          *  Cpi(k)                 / (latent_heat * dt)  ! original method: taking the lower layer value
                        bottom_melt = bottom_melt + Delta_T * (zeta(k) - zeta(k-1)) * Hi * 0.50 * (Cpi(k) + Cpi(k-1)) / (latent_heat * dt)  ! new      method: taking the average value of the zeta layer
                    END IF

                    ! Debug check on extreme values:
                    !IF(bottom_melt < -1.0) write(*,*) i,j, k, bottom_melt, Delta_T, Hi, Ti_new(k), Ti_pmp(k), Ti_new(k-1), Ti_pmp(k-1), Ti_new(k) - Ti_pmp(k), Ti_new(k-1) - Ti_pmp(k-1)

                    ! Adjust those 3D grid points on internal ice layers which still have a temperature above Ti_pmp (see equation (X.45)):
                    Ti_new(k) = Ti_pmp(k)
                !ELSE IF(Ti_new(k) > Ti_pmp(k)) THEN 
                ! If this ELSE IF-block is not in use, it means that the bottom Ti(k-1) (of a layer k-1/2) is only adjusted to Ti_pmp if on average (at k-1/2) this layer is above Ti_pmp(k-1/2)
                ! ! Checking for the case that the layer average Delta_T < 0 but still the lowest layer has Ti_new > Ti_pmp.
                ! ! In that case here the temperature is set to the pressure melting point temperature but no melt is added.
                ! Ti_new(k) = Ti_pmp(k)
                ELSE IF(Ti_new(k) < registration_limit_low_Ti) THEN
                    WRITE(*,'(3(A, I4), 2(A, F12.4), A)') &
                    ' WARNING: A very low temperature is detected: Ti_new(', k, ',', j, ',', i, ') = ', Ti_new(k), '. At this point the vertical temperature profile is replaced by the Robin solution.'
                    IF(.FALSE.) THEN
                        WRITE(*,'(A)') '    i    j  is_float   Ts           dHs_dx      dHs_dy       Hi        zeta_z     Ti_pmp reached'
                        WRITE(*,'( 3I4, 5F12.5)') i, j, is_float, Ts, dHs_dx, dHs_dy, Hi, zeta_z
                        WRITE(*,'(A)') '   k     Ti_pmp      Ti_new        Ti       Cpi         Ki            zeta_t       zeta_x       zeta_y'
                        WRITE(*,'(I4, 4F12.4, F13.1, F12.4, 2F14.9)') k, Ti_pmp(k), Ti_new(k), Ti(k), Cpi(k), Ki(k), zeta_t(k), zeta_x(k), zeta_y(k)
                        WRITE(*,'(A)') ''
                    END IF

                    !Ti_new(:) = Ts + zeta(:) * (seawater_temperature - Ts)  ! Alternative at shelf: Using the seawater temperature
                    !Ti_new(:) = Ti(:)                                           ! Alternative at sheet: Using the vertical profile of the previous time step for this point
                    Ti_new(:) = robin_solution(zeta,Ts,Hi,mb_net,q_bottom,is_float,rho_ice)
                    
                END IF
            END DO

        END IF

        bottom_melt = MAX(bottom_melt, maximum_bottom_melt)

        ! Replace solution with robin_solution to test it 
        !Ti_new(:) = robin_solution(zeta,Ts,Hi,mb_net,q_bottom,is_float,rho_ice)
        !bottom_melt = 0.0 

        return 
    
    end subroutine temperature_imau

    subroutine calc_zeta_helpers(help,zeta,dx,dy)

        implicit none 

        type(zeta_helper_type), intent(INOUT) :: help 
        real(prec), intent(IN)                :: zeta(:) 
        real(prec), intent(IN)                :: dx 
        real(prec), intent(IN)                :: dy

        ! Local variables 
        integer     :: nz, k  
        real(prec)  :: a, b, c

        nz = size(zeta,1)

        if (allocated(help%a_zeta)) then 
            ! Deallocate all object values 
            deallocate(help%a_zeta)
            deallocate(help%a_zetazeta)
            deallocate(help%a_zeta_minus)
            deallocate(help%b_zeta)
            deallocate(help%b_zetazeta)
            deallocate(help%b_zeta_minus)
            deallocate(help%c_zeta)
            deallocate(help%c_zetazeta)
            deallocate(help%z_zeta_minus)

        end if 

        ! Allocate to match length of zeta 
        allocate(help%a_zeta(nz))
        allocate(help%a_zetazeta(nz))
        allocate(help%a_zeta_minus(nz))
        allocate(help%b_zeta(nz))
        allocate(help%b_zetazeta(nz))
        allocate(help%b_zeta_minus(nz))
        allocate(help%c_zeta(nz))
        allocate(help%c_zetazeta(nz))
        allocate(help%z_zeta_minus(nz))

        help%a_zeta         = 0.0 
        help%a_zetazeta     = 0.0 
        help%a_zeta_minus   = 0.0 
        help%b_zeta         = 0.0 
        help%b_zetazeta     = 0.0 
        help%b_zeta_minus   = 0.0 
        help%c_zeta         = 0.0 
        help%c_zetazeta     = 0.0  
        help%z_zeta_minus   = 0.0 

        do k = 2, nz-1 

            a = zeta(k)   - zeta(k-1)
            b = zeta(k+1) - zeta(k) 

            help%a_zeta(k)         =   -b  / (a*(a+b))
            help%a_zetazeta(k)     =  2.0  / (a*(a+b))
            help%b_zeta(k)         = (b-a) / (a*b)
            help%b_zetazeta(k)     = -2.0  /  (a*b)
            help%c_zeta(k)         =    a  / (b*(a+b))
            help%c_zetazeta(k)     =  2.0  / (b*(a+b))
            
        end do  

        k = nz 
        a = zeta(k)   - zeta(k-1) 
        c = zeta(k)   - zeta(k-2)
        help%a_zeta_minus(k)   =   c   / (a*(a-c))
        help%b_zeta_minus(k)   = (a+c) / (a*c)
        help%z_zeta_minus(k)   =   a   / (c*(c-a))

        return 

    end subroutine calc_zeta_helpers 

    SUBROUTINE calculate_zeta_derivatives(t,x,y,z, Hi, dHi_dt, dHi_dx, dHi_dy, dHs_dt, dHs_dx, dHs_dy,zeta)
        ! This subroutine calculates the derivatives of zeta, which are used in the transformation to the t,x,y,zeta-coordinates. 
        ! This zeta derivatives are the Jacobians. See table 3.
        
        implicit none 

        real(prec), intent(OUT) :: t(:) 
        real(prec), intent(OUT) :: x(:) 
        real(prec), intent(OUT) :: y(:) 
        real(prec), intent(OUT) :: z 
        
        real(prec), intent(IN)  :: Hi
        real(prec), intent(IN)  :: dHi_dt
        real(prec), intent(IN)  :: dHi_dx
        real(prec), intent(IN)  :: dHi_dy
        real(prec), intent(IN)  :: dHs_dt
        real(prec), intent(IN)  :: dHs_dx
        real(prec), intent(IN)  :: dHs_dy
        real(prec), intent(IN)  :: zeta(:) 

        ! Local variables:
        integer    :: k, nz 
        real(prec) :: inverse_Hi                ! Contains the inverse of Hi

        nz = size(t,1) 

        inverse_Hi = 1.0 / Hi

        do k = 1, nz
            t(k) = inverse_Hi * (dHs_dt  - zeta(k) * dHi_dt)
            x(k) = inverse_Hi * (dHs_dx  - zeta(k) * dHi_dx)
            y(k) = inverse_Hi * (dHs_dy  - zeta(k) * dHi_dy)
        end do

        z  = -inverse_Hi
            
        return 

    end subroutine calculate_zeta_derivatives

    function robin_solution(zeta, Ts, Hi, mb_net, ghf, is_float, rho_ice) result(Ti)
        ! This function calculates for one horizontal grid point the temperature profiles
        ! using the surface temperature and the geothermal heat flux as boundary conditions.
        ! See Robin solution in: Cuffey & Paterson 2010, 4th ed, chapter 9, eq. (9.13) - (9.22).
        
        implicit none 

        ! Input variables:
        real(prec), intent(IN)  :: zeta(:)              ! Vertical scale (0-1)
        real(prec), intent(IN)  :: Ts                   ! Surface temperature [K]
        real(prec), intent(IN)  :: Hi                   ! Ice thickness [m]
        real(prec), intent(IN)  :: mb_net               ! Netto mass balance [meter ice equivalent per year]
        real(prec), intent(IN)  :: ghf                  ! [J a-1 m-2] Geothermal heat flux 
        logical,    intent(IN)  :: is_float             ! is point floating (true) or grounded (false)
        real(prec), intent(IN)  :: rho_ice             ! [kg m-3] Ice density 
        
        ! Result variables:
        real(prec) :: Ti(size(zeta,1))                ! The calculated new ice temperature distribution

        ! Local variables:
        integer :: k, nz 
        real(prec)    thermal_length_scale                               ! thermal length scale [m] (the z*)
        real(prec)    distance_above_bed                                 ! distance above bed [m]
        real(prec)    error_function_1                                   ! error function
        real(prec)    error_function_2                                   ! error function

        real(prec), parameter :: Hi_min                = 0.1           ! [m]
        real(prec), parameter :: Claus_Clap_gradient   = 8.7E-04       ! [K/m]
        real(prec), parameter :: triple_point_of_water = 273.16        ! [K]
        real(prec), parameter :: delta_Ti_pmp_bottom   = 0.0           ! [K]
        real(prec), parameter :: seawater_temperature  = 271.15        ! [K]

        real(prec), parameter :: c_0_specific_heat        = 2127.5
        real(prec), parameter :: kappa_0_ice_conductivity =  9.828
        real(prec), parameter :: kappa_e_ice_conductivity =  5.7e-3
        real(prec), parameter :: thermal_conductivity_robin  = kappa_0_ice_conductivity * exp(-kappa_e_ice_conductivity * triple_point_of_water)
        

        real(prec) :: thermal_diffusivity_robin
        real(prec) :: bottom_temperature_gradient_robin

        thermal_diffusivity_robin = thermal_conductivity_robin / (rho_ice * c_0_specific_heat)
        bottom_temperature_gradient_robin = -ghf / thermal_conductivity_robin

        nz = size(Ti,1) 

        ! Calculation of temperature profile:
        IF(Hi > Hi_min .AND. mb_net > 0 .AND. (.not. is_float)) THEN
            ! The Robin solution can be used to estimate the subsurface temperature profile
            thermal_length_scale = sqrt(2.0 * thermal_diffusivity_robin * Hi / mb_net)
            DO k = 1, nz
                distance_above_bed = (1.0 - zeta(k)) * Hi
                error_function_1 = error_function(distance_above_bed / thermal_length_scale)
                error_function_2 = error_function(Hi / thermal_length_scale)
                Ti(k) = Ts + sqrt(pi) / 2.0 * thermal_length_scale * bottom_temperature_gradient_robin * (error_function_1 - error_function_2)
            END DO
        ELSE IF(mb_net <= 0.0 .AND. Hi > Hi_min) THEN
            ! Ablation area: use linear temperature profile from Ts to (offset below) T_pmp
            Ti(:) = Ts + ((triple_point_of_water - Claus_Clap_gradient * Hi - delta_Ti_pmp_bottom) - Ts) * zeta(:)
        ELSE IF(is_float .AND. Hi > Hi_min) THEN
            ! This option initializes the temperature at shelf points linear between Ts and seawater temperature:
            Ti(:) = Ts + zeta(:) * (seawater_temperature - Ts)
        ELSE
            ! Only very thin layer: use Ts for entire ice profile
            Ti(:) = Ts
        END IF

        ! Correct all temperatures above T_pmp:
        do k = 1, nz
            IF(Ti(k) > triple_point_of_water - Claus_Clap_gradient * Hi * zeta(k)) &
                Ti(k) = triple_point_of_water - Claus_Clap_gradient * Hi * zeta(k)
        END do

        return
    
    end function robin_solution

end module icetemp_imau
