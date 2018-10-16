module icetemp_mali 
    ! Wrapping the MALI icetemp solution for yelmo 

    use defs, only : prec, pi, g, sec_year, T0, rho_ice, rho_sw, rho_w 
    use solver_tridiagonal, only : solve_tridiag 
    use thermodynamics, only : calc_advec_horizontal_column_sico1, calc_advec_horizontal_column, &
                                calc_T_pmp

    implicit none
    
    ! Impose parameter choice here for testing conductive bedrock or not
    logical,    parameter :: conductive_bedrock = .FALSE. 
    
    private
    public :: calc_icetemp_mali_3D_up

    public :: calc_mali_temp_column
    public :: calc_sigt_terms



contains 

    subroutine calc_icetemp_mali_3D_up(T_ice,T_pmp,cp,ct,ux,uy,uz,Q_strn,Q_b, &
                                            Q_geo,T_srf,H_ice,H_w,smb,bmb,f_grnd,sigma,sigt,dsigt_a,dsigt_b,dt,dx)
        ! Solver for thermodynamics of ice 
        ! Note sigma=height, k=1 base, k=nz surface 
        ! Note T_ice, T_pmp in [degC], not [K]

        implicit none 

        real(prec), intent(OUT) :: T_ice(:,:,:)     ! [degC] Ice column temperature
        real(prec), intent(IN)  :: T_pmp(:,:,:)     ! [degC] Pressure melting point temp.
        real(prec), intent(IN)  :: cp(:,:,:)        ! [J kg-1 K-1] Specific heat capacity
        real(prec), intent(IN)  :: ct(:,:,:)        ! [J a-1 m-1 K-1] Heat conductivity 
        real(prec), intent(IN)  :: ux(:,:,:)        ! [m a-1] Horizontal x-velocity 
        real(prec), intent(IN)  :: uy(:,:,:)        ! [m a-1] Horizontal y-velocity 
        real(prec), intent(IN)  :: uz(:,:,:)        ! [m a-1] Vertical velocity 
        real(prec), intent(IN)  :: Q_strn(:,:,:)    ! [K a-1] Internal strain heat production in ice
        real(prec), intent(IN)  :: Q_b(:,:)         ! [J a-1 m-2] Basal frictional heat production 
        real(prec), intent(IN)  :: Q_geo(:,:)       ! [mW m-2] Geothermal heat flux 
        real(prec), intent(IN)  :: T_srf(:,:)       ! [degC] Surface temperature 
        real(prec), intent(IN)  :: H_ice(:,:)       ! [m] Ice thickness 
        real(prec), intent(IN)  :: H_w(:,:)         ! [m] Basal water layer thickness 
        real(prec), intent(IN)  :: smb(:,:)         ! [m a-1] Surface mass balance (melting is negative)
        real(prec), intent(IN)  :: bmb(:,:)         ! [m a-1] Basal mass balance (melting is negative)
        real(prec), intent(IN)  :: f_grnd(:,:)      ! [--] Floating point or grounded?
        real(prec), intent(IN)  :: sigma(:)         ! [--] Vertical sigma coordinates (sigma==height), dynamics
        real(prec), intent(IN)  :: sigt(:)          ! [--] Vertical sigma coordinates (sigma==height), thermodynamics
        real(prec), intent(IN)  :: dsigt_a(:)       ! d Vertical height axis (0:1) 
        real(prec), intent(IN)  :: dsigt_b(:)       ! d Vertical height axis (0:1) 
        real(prec), intent(IN)  :: dt               ! [a] Time step 
        real(prec), intent(IN)  :: dx               ! [a] Horizontal grid step 
        
        ! Local variables
        integer :: i, j, k, nx, ny, nz, nzt   
        real(prec), allocatable  :: advecxy(:)   ! [K a-1 m-2] Horizontal heat advection 
        logical :: is_float 
        real(prec) :: H_w_dot 
        real(prec) :: ghf_conv, Q_geo_now 

        ghf_conv = 1e-3*sec_year   ! [mW m-2] => [J m-2 a-1]

        nx  = size(T_ice,1)
        ny  = size(T_ice,2)
        nzt = size(T_ice,3)
        nz  = size(uz,3)

        allocate(advecxy(nzt))
        advecxy = 0.0 

        do j = 3, ny-2
        do i = 3, nx-2 
        
            ! Get geothermal heat flux in proper units 
            Q_geo_now = Q_geo(i,j)*ghf_conv 

            ! Determine if point is floating 
            is_float = (f_grnd(i,j) .eq. 0.0)
            
            if (H_ice(i,j) .gt. 2.0) then 
                ! Thick ice exists, call thermodynamic solver for the column

                ! Pre-calculate the contribution of horizontal advection to column solution
                call calc_advec_horizontal_column(advecxy,T_ice,ux,uy,dx,i,j)
                
                call calc_mali_temp_column(T_ice(i,j,:),T_pmp(i,j,:),cp(i,j,:),ct(i,j,:),uz(i,j,:), &
                                                Q_strn(i,j,:),advecxy,Q_b(i,j),Q_geo(i,j),T_srf(i,j), &
                                                H_ice(i,j),H_w(i,j),bmb(i,j),is_float, &
                                                sigma,sigt,dsigt_a,dsigt_b,dt)

            else 
                ! Ice is too thin, prescribe ice temperature for now

                ! To do 
                T_ice(i,j,:) = T_pmp(i,j,:) 

            end if 


        end do 
        end do 

        return 

    end subroutine calc_icetemp_mali_3D_up

    subroutine calc_mali_temp_column(T_ice,T_pmp,cp,ct,uz,Q_strn,advecxy,Q_b, &
                                            Q_geo,T_srf,H_ice,H_w,bmb,is_float, &
                                            sigma,sigt,dsigt_a,dsigt_b,dt)
        ! Thermodynamics solver for a given column of ice 
        ! Note sigma=height, k=1 base, k=nz surface 
        ! Note T_ice, T_pmp in [degC], not [K]
        ! Note: nz = number of vertical boundaries (including sigma=0.0 and sigma=1.0), 
        ! temperature is defined for cell centers, plus a value at the surface and the base
        ! so nzt = nz + 1 

        implicit none 

        real(prec), intent(INOUT) :: T_ice(:)     ! nz+1 [degC] Ice column temperature
        real(prec), intent(IN)    :: T_pmp(:)     ! nz+1 [degC] Pressure melting point temp.
        real(prec), intent(IN)    :: cp(:)        ! nz+1 [J kg-1 K-1] Specific heat capacity
        real(prec), intent(IN)    :: ct(:)        ! nz+1 [J a-1 m-1 K-1] Heat conductivity 
        real(prec), intent(IN)    :: uz(:)        ! nz   [m a-1] Vertical velocity 
        real(prec), intent(IN)    :: Q_strn(:)    ! nz+1 [K a-1] Internal strain heat production in ice
        real(prec), intent(IN)    :: advecxy(:)   ! nz+1 [K a-1] Horizontal heat advection 
        real(prec), intent(IN)    :: Q_b          ! [J a-1 m-2] Basal frictional heat production 
        real(prec), intent(IN)    :: Q_geo        ! [mW m-2] Geothermal heat flux 
        real(prec), intent(IN)    :: T_srf        ! [degC] Surface temperature 
        real(prec), intent(IN)    :: H_ice        ! [m] Ice thickness 
        real(prec), intent(IN)    :: H_w          ! [m] Basal water layer thickness 
        real(prec), intent(IN)    :: bmb          ! [m a-1] Basal mass balance (melting is negative)
        logical,    intent(IN)    :: is_float     ! [--] Floating point or grounded?
        real(prec), intent(IN)    :: sigma(:)     ! nz   [--] Vertical sigma coordinates (sigma==height)
        real(prec), intent(IN)    :: sigt(:)      ! nz+1 [--] Vertical height axis temperature (0:1) 
        real(prec), intent(IN)    :: dsigt_a(:)   ! nz+1 [--] d Vertical height axis (0:1) 
        real(prec), intent(IN)    :: dsigt_b(:)   ! nz+1 [--] d Vertical height axis (0:1) 
        real(prec), intent(IN)    :: dt           ! [a] Time step 

        ! Local variables 
        integer :: k, nz, nzt, ki  
        real(prec) :: T_pmp_base
        real(prec) :: Q_geo_now, ghf_conv  
        !real(prec), allocatable :: advecz_aa(:) ! nzt 

        real(prec), allocatable :: subd(:)     ! nzt 
        real(prec), allocatable :: diag(:)     ! nzt  
        real(prec), allocatable :: supd(:)     ! nzt 
        real(prec), allocatable :: rhs(:)      ! nzt 
        real(prec), allocatable :: solution(:) ! nzt
        real(prec) :: T_base, dsigmaBot, fac 
        
        nz  = size(sigma,1)
        nzt = size(sigt,1)   ! == nz+1, base point, layer midpoints where T is defined, and surface point 

        allocate(subd(nzt))
        allocate(diag(nzt))
        allocate(supd(nzt))
        allocate(rhs(nzt))
        allocate(solution(nzt))

        ! Get pressure melting point temperature at the base [celcius]
        !T_pmp_base = calc_T_pmp(H_ice,sigma=0.0_prec,T0=0.0_prec)
        T_pmp_base = T_pmp(1) 

        ! Copy basal temperature for convenience 
        T_base     = T_ice(1) 

        ! Get geothermal heat flux in proper units 
        Q_geo_now = Q_geo*1e-3*sec_year   ! [mW m-2] => [J m-2 a-1]

        ! Step 1: apply vertical advection (for explicit testing)
        !allocate(advecz_aa(nzt))
        !advecz_aa = 0.0
        !call advection_1D_upwind(advecz_aa,T_ice,uz,H_ice,sigt)
        !T_ice = T_ice - dt*advecz_aa 

        ! Step 2: apply vertical implicit diffusion-advection 
        
        ! Ice base
        if (is_float) then
            ! Floating ice - set temperature equal to basal temperature 

            subd(1) = 0.0_prec
            diag(1) = 1.0_prec
            supd(1) = 0.0_prec
            rhs(1)  = T_base 

        else 
            ! Grounded ice 

!             if (abs(T_base - T_pmp_base) < 0.001_prec) then
            if (T_base > (T_pmp_base - 0.001_prec)) then
                ! Temperate at bed 
                ! Hold basal temperature at pressure melting point

                subd(1) = 0.0_prec
                diag(1) = 1.0_prec
                supd(1) = 0.0_prec
                rhs(1)  = T_pmp_base 

            else   
                ! Frozen at bed
                ! maintain balance of heat sources and sinks
                ! (conductive flux, geothermal flux, and basal friction)

                ! Note: basalHeatFlux is generally >= 0, since defined as positive up

                ! calculate dsigma for the bottom layer between the basal boundary and the temperature point above
                dsigmaBot = sigt(2) - sigt(1) 

                ! backward Euler flux basal boundary condition
                subd(1) =  0.0_prec
                diag(1) =  1.0_prec
                supd(1) = -1.0_prec
                rhs(1)  = (Q_b + Q_geo_now) * dsigmaBot*H_ice / ct(1)

            end if   ! melting or frozen

        end if  ! floating or grounded 

        ! Ice interior, layers 1:nzt  (matrix elements 2:nzt+1)

        do k = 2, nzt-1

            ! No advection (diffusion only)
!             fac     = dt * ct(k) / (rho_ice*cp(k)) / H_ice**2
!             subd(k) = -fac * dsigt_a(k)
!             supd(k) = -fac * dsigt_b(k)
!             diag(k) = 1.0_prec - (subd(k) + supd(k))
!             rhs(k)  = T_ice(k) + dt*Q_strn(k) 
            
            ! With implicit advection (diffusion + advection)
            fac     = dt * ct(k) / (rho_ice*cp(k)) / H_ice**2
            subd(k) = -fac*dsigt_a(k) - 0.5*(uz(k-1)+uz(k))*dt / (2.0*H_ice*(sigma(k)-sigma(k-1)))
            supd(k) = -fac*dsigt_b(k) + 0.5*(uz(k-1)+uz(k))*dt / (2.0*H_ice*(sigma(k)-sigma(k-1)))
            diag(k) = 1.0_prec - (-fac*dsigt_a(k)) - (-fac*dsigt_b(k))
            rhs(k)  = T_ice(k) + dt*Q_strn(k) - dt*advecxy(k) 

        end do 

        ! Ice surface 
        subd(nzt) = 0.0_prec
        diag(nzt) = 1.0_prec
        supd(nzt) = 0.0_prec
        rhs(nzt)  = T_srf

        ! Call solver 
        call solve_tridiag(subd,diag,supd,rhs,solution)

        ! Copy the solution into the temperature variables
        T_base = solution(1)
        T_ice  = solution
        !T_srf = solution(nzt)
        
        return 

    end subroutine calc_mali_temp_column

    subroutine advection_1D_upwind(advecz,Q,uz,H_ice,sigt)
        ! Calculate vertical advection term advecz, which enters
        ! advection equation as
        ! Q_new = Q - dt*advecz = Q - dt*u*dQ/dx

        implicit none 

        real(prec), intent(OUT)   :: advecz(:) ! nzt: bottom, cell centers, top 
        real(prec), intent(INOUT) :: Q(:)      ! nzt: bottom, cell centers, top 
        real(prec), intent(IN)    :: uz(:)     ! nzt-1 == nz, cell boundaries
        real(prec), intent(IN)    :: H_ice     ! Ice thickness 
        real(prec), intent(IN)    :: sigt(:)   ! nzt, cell centers
        
        ! Local variables
        integer :: k, nzt   
        real(prec) :: u_aa, dx  

        nzt = size(sigt,1)

        advecz = 0.0 

        ! Loop over internal cell centers and perform upwind advection 
        do k = 2, nzt-1 
            u_aa = 0.5*(uz(k-1)+uz(k)) ! Get velocity at cell center
            if (u_aa < 0.0) then 
                ! Upwind negative
                dx = H_ice*(sigt(k+1)-sigt(k))
                advecz(k) = uz(k)*(Q(k+1)-Q(k))/dx   
            else
                ! Upwind positive
                dx = H_ice*(sigt(k)-sigt(k-1))
                advecz(k) = uz(k-1)*(Q(k)-Q(k-1))/dx
            end if 
        end do 

        return 

    end subroutine advection_1D_upwind

    subroutine calc_sigt_terms(sigt,dsig_a,dsig_b,sigma)
        ! sigma = depth axis (1: base, nz: surface)
        ! Calculate ak, bk terms as defined in Hoffmann et al (2018)
        implicit none 

        real(prec), intent(INOUT) :: sigt(:)      ! nz+1 
        real(prec), intent(INOUT) :: dsig_a(:)    ! nz+1
        real(prec), intent(INOUT) :: dsig_b(:)    ! nz+1
        real(prec), intent(IN)    :: sigma(:)     ! nz 

        ! Local variables 
        integer :: k, nz_layers, nz, nzt   

        nz  = size(sigma)
        nzt = size(sigt)

        ! Get sigt (midpoints of sigma layers - between sigma values)
        ! Note: that first layer and last layer are not equally spaced with the remaining vector
        sigt(1) = 0.0 
        do k = 1, nzt-2
            sigt(k+1) = 0.5 * (sigma(k+1) + sigma(k))
        end do 
        sigt(nzt) = 1.0 

        ! Initialize dsig_a/dsig_b to zero, first and last indices will not be used (end points)
        dsig_a = 0.0 
        dsig_b = 0.0 

        do k = 2, nzt-2 
            dsig_a(k) = 1.0/ ( (sigma(k) - sigma(k-1)) * &
                                (sigt(k) - sigt(k-1)) )
        enddo

        do k = 2, nzt-2
            dsig_b(k) = 1.0/ ( (sigma(k) - sigma(k-1)) * &
                                (sigt(k+1) - sigt(k)) )
        end do

        return 

    end subroutine calc_sigt_terms

end module icetemp_mali


