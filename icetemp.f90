module icetemp 
    ! Module contains the ice temperature and basal mass balance (grounded) solution

    use defs, only : prec, pi, g, sec_year, T0, rho_ice, rho_sw, rho_w, L_ice  
    use solver_tridiagonal, only : solve_tridiag 
    use thermodynamics, only : calc_bmb_grounded, calc_advec_vertical_column, calc_advec_horizontal_column, &
                                calc_T_pmp, calc_T_base_shlf_approx, calc_temp_linear_column

    implicit none
    
    private
    public :: calc_icetemp_3D
    public :: calc_temp_column 
    public :: calc_dzeta_terms

contains 

    subroutine calc_icetemp_3D(T_ice,bmb_grnd,dTdz_b,T_pmp,cp,ct,ux,uy,uz,Q_strn,Q_b,Q_geo, &
                            T_srf,H_ice,H_w,f_grnd,zeta_aa,zeta_ac,dzeta_a,dzeta_b,dt,dx)
        ! Solver for thermodynamics of ice 
        ! Note zeta=height, k=1 base, k=nz surface 

        implicit none 

        real(prec), intent(INOUT) :: T_ice(:,:,:)   ! [K] Ice column temperature
        real(prec), intent(INOUT) :: bmb_grnd(:,:)  ! [m a-1] Basal mass balance (melting is negative)
        real(prec), intent(OUT)   :: dTdz_b(:,:)    ! [K m-1] Ice temperature gradient at the base
        real(prec), intent(IN)    :: T_pmp(:,:,:)   ! [K] Pressure melting point temp.
        real(prec), intent(IN)    :: cp(:,:,:)      ! [J kg-1 K-1] Specific heat capacity
        real(prec), intent(IN)    :: ct(:,:,:)      ! [J a-1 m-1 K-1] Heat conductivity 
        real(prec), intent(IN)    :: ux(:,:,:)      ! [m a-1] Horizontal x-velocity 
        real(prec), intent(IN)    :: uy(:,:,:)      ! [m a-1] Horizontal y-velocity 
        real(prec), intent(IN)    :: uz(:,:,:)      ! [m a-1] Vertical velocity 
        real(prec), intent(IN)    :: Q_strn(:,:,:)  ! [K a-1] Internal strain heat production in ice
        real(prec), intent(IN)    :: Q_b(:,:)       ! [J a-1 m-2] Basal frictional heat production 
        real(prec), intent(IN)    :: Q_geo(:,:)     ! [mW m-2] Geothermal heat flux 
        real(prec), intent(IN)    :: T_srf(:,:)     ! [K] Surface temperature 
        real(prec), intent(IN)    :: H_ice(:,:)     ! [m] Ice thickness 
        real(prec), intent(IN)    :: H_w(:,:)       ! [m] Basal water layer thickness 
        real(prec), intent(IN)    :: f_grnd(:,:)    ! [--] Floating point or grounded?
        real(prec), intent(IN)    :: zeta_aa(:)     ! [--] Vertical sigma coordinates (zeta==height), aa-nodes
        real(prec), intent(IN)    :: zeta_ac(:)     ! [--] Vertical sigma coordinates (zeta==height), ac-nodes
        real(prec), intent(IN)    :: dzeta_a(:)    ! d Vertical height axis (0:1) 
        real(prec), intent(IN)    :: dzeta_b(:)    ! d Vertical height axis (0:1) 
        real(prec), intent(IN)    :: dt             ! [a] Time step 
        real(prec), intent(IN)    :: dx             ! [a] Horizontal grid step 
        
        ! Local variables
        integer :: i, j, k, nx, ny, nz_aa, nz_ac  
        real(prec), allocatable  :: advecxy(:)   ! [K a-1 m-2] Horizontal heat advection 
        logical :: is_float 
        real(prec) :: H_w_dot 

        nx    = size(T_ice,1)
        ny    = size(T_ice,2)
        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1)

        allocate(advecxy(nz_aa))
        advecxy = 0.0 

        do j = 3, ny-2
        do i = 3, nx-2 
            
            ! Determine if point is floating 
            is_float = (f_grnd(i,j) .eq. 0.0)
            
            if (H_ice(i,j) .gt. 10.0) then 
                ! Thick ice exists, call thermodynamic solver for the column

                ! Pre-calculate the contribution of horizontal advection to column solution
                call calc_advec_horizontal_column(advecxy,T_ice,ux,uy,dx,i,j)
                
                call calc_temp_column(T_ice(i,j,:),bmb_grnd(i,j),dTdz_b(i,j),T_pmp(i,j,:),cp(i,j,:),ct(i,j,:), &
                                        uz(i,j,:),Q_strn(i,j,:),advecxy,Q_b(i,j),Q_geo(i,j),T_srf(i,j),H_ice(i,j), &
                                        H_w(i,j),is_float,zeta_aa,zeta_ac,dzeta_a,dzeta_b,dt)

            else ! H_ice(i,j) .le. 10.0
                ! Ice is too thin or zero, prescribe linear temperature profile
                ! between temperate ice at base and surface temperature 

                T_ice(i,j,:) = calc_temp_linear_column(T_srf(i,j),T_pmp(i,j,1),T_pmp(i,j,nz_aa),zeta_aa)

            end if 

        end do 
        end do 

        ! Fill in borders 
        T_ice(2,:,:)    = T_ice(3,:,:) 
        T_ice(1,:,:)    = T_ice(3,:,:) 
        T_ice(nx-1,:,:) = T_ice(nx-2,:,:) 
        T_ice(nx,:,:)   = T_ice(nx-2,:,:) 
        
        T_ice(:,2,:)    = T_ice(:,3,:) 
        T_ice(:,1,:)    = T_ice(:,3,:) 
        T_ice(:,ny-1,:) = T_ice(:,ny-2,:) 
        T_ice(:,ny,:)   = T_ice(:,ny-2,:) 
        
        ! Without conductive bedrock: fill in bedrock solution with steady-state diffusion solution
        ! (linear profile with slope of Q_geo/k)

        ! TO DO 

!         do j = 1, ny
!         do i = 1, nx 
            
!             ! Prescribe solution for bedrock points:
!             ! Calculate temperature in the bedrock as a linear gradient of ghf
!             ! Note: Using T(nzb) from the previous timestep for simplicity
!             do k = 1, nzm
!                 T_rock(i,j,k) = T(nzb)+dzm*(nzm-k+1)*Q_geo_now/cm
!             end do
        
!         end do 
!         end do 

        return 

    end subroutine calc_icetemp_3D

    subroutine calc_temp_column(T_ice,bmb_grnd,dTdz_b,T_pmp,cp,ct,uz,Q_strn,advecxy,Q_b,Q_geo, &
                                T_srf,H_ice,H_w,is_float,zeta_aa,zeta_ac,dzeta_a,dzeta_b,dt)
        ! Thermodynamics solver for a given column of ice 
        ! Note zeta=height, k=1 base, k=nz surface 
        ! Note: nz = number of vertical boundaries (including zeta=0.0 and zeta=1.0), 
        ! temperature is defined for cell centers, plus a value at the surface and the base
        ! so nz_ac = nz_aa - 1 

        ! For notes on implicit form of advection terms, see eg http://farside.ph.utexas.edu/teaching/329/lectures/node90.html
        
        implicit none 

        real(prec), intent(INOUT) :: T_ice(:)     ! nz_aa [K] Ice column temperature
        real(prec), intent(INOUT) :: bmb_grnd     ! [m a-1] Basal mass balance (melting is negative)
        real(prec), intent(OUT)   :: dTdz_b       ! [K m-1] Basal temperature gradient
        real(prec), intent(IN)    :: T_pmp(:)     ! nz_aa [K] Pressure melting point temp.
        real(prec), intent(IN)    :: cp(:)        ! nz_aa [J kg-1 K-1] Specific heat capacity
        real(prec), intent(IN)    :: ct(:)        ! nz_aa [J a-1 m-1 K-1] Heat conductivity 
        real(prec), intent(IN)    :: uz(:)        ! nz_ac [m a-1] Vertical velocity 
        real(prec), intent(IN)    :: Q_strn(:)    ! nz_aa [K a-1] Internal strain heat production in ice
        real(prec), intent(IN)    :: advecxy(:)   ! nz_aa [K a-1] Horizontal heat advection 
        real(prec), intent(IN)    :: Q_b          ! [J a-1 m-2] Basal frictional heat production 
        real(prec), intent(IN)    :: Q_geo        ! [mW m-2] Geothermal heat flux 
        real(prec), intent(IN)    :: T_srf        ! [K] Surface temperature 
        real(prec), intent(IN)    :: H_ice        ! [m] Ice thickness 
        real(prec), intent(IN)    :: H_w          ! [m] Basal water layer thickness 
        logical,    intent(IN)    :: is_float     ! [--] Floating point or grounded?
        real(prec), intent(IN)    :: zeta_aa(:)  ! nz_aa [--] Vertical sigma coordinates (zeta==height), layer centered aa-nodes
        real(prec), intent(IN)    :: zeta_ac(:)  ! nz_ac [--] Vertical height axis temperature (0:1), layer edges ac-nodes
        real(prec), intent(IN)    :: dzeta_a(:)  ! nz_aa [--] Solver discretization helper variable ak
        real(prec), intent(IN)    :: dzeta_b(:)  ! nz_aa [--] Solver discretization helper variable bk
        real(prec), intent(IN)    :: dt           ! [a] Time step 

        ! Local variables 
        integer :: k, nz_aa, nz_ac
        real(prec) :: Q_geo_now, ghf_conv 
        real(prec) :: H_w_dot
        real(prec) :: melt_internal, T_excess  

        real(prec), allocatable :: advecz(:)   ! nz_aa, for explicit vertical advection solving
        logical, parameter      :: test_expl_advecz = .FALSE. 

        real(prec), allocatable :: subd(:)     ! nz_aa 
        real(prec), allocatable :: diag(:)     ! nz_aa  
        real(prec), allocatable :: supd(:)     ! nz_aa 
        real(prec), allocatable :: rhs(:)      ! nz_aa 
        real(prec), allocatable :: solution(:) ! nz_aa
        real(prec) :: dzetaBot, fac, fac_a, fac_b, uz_aa, dz  

        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1)

        allocate(subd(nz_aa))
        allocate(diag(nz_aa))
        allocate(supd(nz_aa))
        allocate(rhs(nz_aa))
        allocate(solution(nz_aa))

        ! Get geothermal heat flux in proper units 
        Q_geo_now = Q_geo*1e-3*sec_year   ! [mW m-2] => [J m-2 a-1]

        ! Step 1: apply vertical advection (for explicit testing)
        if (test_expl_advecz) then 
            allocate(advecz(nz_aa))
            advecz = 0.0
            call calc_advec_vertical_column(advecz,T_ice,uz,H_ice,zeta_aa)
            T_ice = T_ice - dt*advecz 
        end if 

        ! Step 2: apply vertical implicit diffusion-advection 
        
        ! Ice base
        if (is_float) then
            ! Floating ice - set temperature equal to basal temperature at pressure melting point, or marine freezing temp

            subd(1) = 0.0_prec
            diag(1) = 1.0_prec
            supd(1) = 0.0_prec
            rhs(1)  = min(calc_T_base_shlf_approx(H_ice),T_pmp(1))

        else 
            ! Grounded ice 

            ! Determine expected change in basal water total for this timestep [m]
            H_w_dot = -(bmb_grnd*(rho_w/rho_ice))*dt 


            if (T_ice(1) .lt. T_pmp(1) .or. H_w+H_w_dot .lt. 0.0_prec) then   
                ! Frozen at bed, or about to become frozen 

                ! maintain balance of heat sources and sinks
                ! (conductive flux, geothermal flux, and basal friction)

                ! Note: basalHeatFlux is generally >= 0, since defined as positive up

                ! calculate dzeta for the bottom layer between the basal boundary and the temperature point above
                dzetaBot = zeta_aa(2) - zeta_aa(1) 

                ! backward Euler flux basal boundary condition
                subd(1) =  0.0_prec
                diag(1) =  1.0_prec
                supd(1) = -1.0_prec
                rhs(1)  = (Q_b + Q_geo_now) * dzetaBot*H_ice / ct(1)

            else 
                ! Temperate at bed 
                ! Hold basal temperature at pressure melting point

                subd(1) = 0.0_prec
                diag(1) = 1.0_prec
                supd(1) = 0.0_prec
                rhs(1)  = T_pmp(1) 

            end if   ! melting or frozen

        end if  ! floating or grounded 

        ! Ice interior layers 2:nz_aa-1
        do k = 2, nz_aa-1

            if (test_expl_advecz) then 
                ! No implicit vertical advection (diffusion only)
                uz_aa = 0.0 

            else
                ! With implicit vertical advection (diffusion + advection)
                uz_aa   = 0.5*(uz(k-1)+uz(k))   ! ac => aa nodes
            end if 

            dz      =  H_ice*(zeta_ac(k)-zeta_ac(k-1))
            

            fac     = dt * ct(k) / (rho_ice*cp(k)) / H_ice**2
            fac_a   = -fac*dzeta_a(k)
            fac_b   = -fac*dzeta_b(k)
            subd(k) = fac_a - uz_aa*dt / (2.0*dz)
            supd(k) = fac_b + uz_aa*dt / (2.0*dz)
            diag(k) = 1.0_prec - fac_a - fac_b
            rhs(k)  = T_ice(k) + dt*Q_strn(k) - dt*advecxy(k) 

        end do 

        ! Ice surface 
        subd(nz_aa) = 0.0_prec
        diag(nz_aa) = 1.0_prec
        supd(nz_aa) = 0.0_prec
        rhs(nz_aa)  = min(T_srf,T0)

        ! Call solver 
        call solve_tridiag(subd,diag,supd,rhs,solution)

        ! Copy the solution into the temperature variables
        T_ice  = solution
        
        ! === Treat basal mass balance and high temperatures ===
        
        ! First calculate internal melt (only allow melting, no accretion)
        
        melt_internal = 0.0 

        do k = nz_aa-1, 2, -1 
            ! Descend from surface to base layer (center of layer)

            ! Store temperature difference with pressure melting point (excess energy)
            T_excess = T_ice(k) - T_pmp(k) 

            ! Calculate basal mass balance as sum of all water produced in column 
            if (T_excess .gt. 0.0) then 
                melt_internal = melt_internal - T_excess * H_ice*(zeta_ac(k)-zeta_ac(k-1))*cp(k) / (L_ice * dt) 
            end if 

            ! Reset temperature to below T_pmp if needed
            if (T_ice(k) .gt. T_pmp(k)) T_ice(k) = T_pmp(k)

        end do 

        ! Make sure base is below pmp too (mass/energy balance handled via bmb_grnd calculation externally)
        k = 1 
        if (T_ice(k) .gt. T_pmp(k)) T_ice(k) = T_pmp(k)


        ! Get temperature gradient at ice base
        if (H_ice .gt. 0.0_prec) then 
            dz = H_ice * (zeta_aa(2)-zeta_aa(1))
            dTdz_b = (T_ice(2) - T_ice(1)) / dz 
        else 
            dTdz_b = 0.0_prec 
        end if 
        
        ! Calculate basal mass balance (valid for grounded ice only)
        call calc_bmb_grounded(bmb_grnd,T_ice(1)-T_pmp(1),dTdz_b,ct(1),rho_ice,Q_b,Q_geo_now,is_float)

        ! Include internal melting in bmb_grnd 
        bmb_grnd = bmb_grnd - melt_internal 

        return 

    end subroutine calc_temp_column

    subroutine calc_dzeta_terms(dzeta_a,dzeta_b,zeta_aa,zeta_ac)
        ! zeta_aa  = depth axis at layer centers (plus base and surface values)
        ! zeta_ac  = depth axis (1: base, nz: surface), at layer boundaries
        ! Calculate ak, bk terms as defined in Hoffmann et al (2018)
        implicit none 

        real(prec), intent(INOUT) :: dzeta_a(:)    ! nz_aa
        real(prec), intent(INOUT) :: dzeta_b(:)    ! nz_aa
        real(prec), intent(IN)    :: zeta_aa(:)    ! nz_aa 
        real(prec), intent(IN)    :: zeta_ac(:)    ! nz_ac 

        ! Local variables 
        integer :: k, nz_layers, nz_aa, nz_ac    

        nz_aa = size(zeta_aa)
        nz_ac = size(zeta_ac)

        ! Note: zeta_aa is calculated outside in the main program 

        ! Initialize dzeta_a/dzeta_b to zero, first and last indices will not be used (end points)
        dzeta_a = 0.0 
        dzeta_b = 0.0 

        do k = 2, nz_aa-1 
            dzeta_a(k) = 1.0/ ( (zeta_ac(k) - zeta_ac(k-1)) * (zeta_aa(k) - zeta_aa(k-1)) )
        enddo

        do k = 2, nz_aa-1
            dzeta_b(k) = 1.0/ ( (zeta_ac(k) - zeta_ac(k-1)) * (zeta_aa(k+1) - zeta_aa(k)) )
        end do

        return 

    end subroutine calc_dzeta_terms
    
end module icetemp


