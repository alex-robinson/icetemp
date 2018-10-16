module icetemp_mali 
    ! Wrapping the MALI icetemp solution for yelmo 

    use defs, only : prec, pi, g, sec_year, T0, rho_ice, rho_sw, rho_w 
    use solver_tridiagonal, only : tridiag 
    use thermodynamics, only : calc_advec_horizontal_column_sico1, calc_advec_horizontal_column, &
                                calc_T_pmp

    implicit none
    
    ! Impose parameter choice here for testing conductive bedrock or not
    logical,    parameter :: conductive_bedrock = .FALSE. 
    
    private
    !public :: calc_icetemp_mali_3D_up
    public :: calc_icetemp_mali_column_up
     
    public :: mali_temp_diffusion_column
    public :: calc_sigt_terms



contains 

    subroutine calc_icetemp_mali_column_up(T_ice,bmb,is_float, H_ice, T_srf, advecxy, ux, uy, uz, dzsdx, dzsdy, dzsrfdt, &
                                      dHicedx, dHicedy, dHicedt, cp, kt, Q_strn, Q_b, T_pmp, mb_net, Q_geo, sigma, &
                                      sigt, dsigt_a, dsigt_b, dx, dt)
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
        real(prec), intent(IN)  :: sigt(:)                               ! Vertical height axis (0:1) 
        real(prec), intent(IN)  :: dsigt_a(:)                               ! d Vertical height axis (0:1) 
        real(prec), intent(IN)  :: dsigt_b(:)                               ! d Vertical height axis (0:1) 
        
        real(prec), intent(IN)  :: dx 
        real(prec), intent(IN)  :: dt 

        ! Local variables
        integer    :: i, j, k, nx, ny, nz, nzt   
        real(prec) :: bottom_melt 
        real(prec) :: Q_geo_now 
        real(prec) :: ghf_conv, ghf_now 

        ghf_conv = 1e-3*sec_year   ! [mW m-2] => [J m-2 a-1]

        nz  = size(sigma,1)
        nzt = size(sigt,1) 

        ! Get geothermal heat flux in proper units 
        Q_geo_now = Q_geo*ghf_conv 


        return 

    end subroutine calc_icetemp_mali_column_up


    subroutine mali_temp_diffusion_column(T_ice,T_base,T_pmp,cp,ct,uz,Q_strn,advecxy,Q_b, &
                                            Q_geo,T_srf,H_ice,H_w,bmb,is_float, &
                                            sigma,sigt,dsigt_a,dsigt_b,dt)
        ! Thermodynamics solver for a given column of ice 
        ! Note sigma=height, k=1 base, k=nz surface 
        ! Note T_ice, T_pmp in [degC], not [K]
        ! Note: nz = number of vertical boundaries (including sigma=0.0 and sigma=1.0), 
        ! temperature is defined for cell centers, plus a value at the surface and the base
        ! so nzt = nz + 1 

        implicit none 

        real(prec), intent(INOUT) :: T_ice(:)   ! nz-1 [degC] Ice column temperature
        real(prec), intent(INOUT) :: T_base     ! [degC] Basal temperature
        real(prec), intent(IN)  :: T_pmp(:)     ! nz-1 [degC] Pressure melting point temp.
        real(prec), intent(IN)  :: cp(:)        ! nz   [J kg-1 K-1] Specific heat capacity
        real(prec), intent(IN)  :: ct(:)        ! nz   [J a-1 m-1 K-1] Heat conductivity 
        real(prec), intent(IN)  :: uz(:)        ! nz   [m a-1] Vertical velocity 
        real(prec), intent(IN)  :: Q_strn(:)    ! nz   [K a-1] Internal strain heat production in ice
        real(prec), intent(IN)  :: advecxy(:)   ! nz   [K a-1] Horizontal heat advection 
        real(prec), intent(IN)  :: Q_b          ! [J a-1 m-2] Basal frictional heat production 
        real(prec), intent(IN)  :: Q_geo        ! [mW m-2] Geothermal heat flux 
        real(prec), intent(IN)  :: T_srf        ! [degC] Surface temperature 
        real(prec), intent(IN)  :: H_ice        ! [m] Ice thickness 
        real(prec), intent(IN)  :: H_w          ! [m] Basal water layer thickness 
        real(prec), intent(IN)  :: bmb          ! [m a-1] Basal mass balance (melting is negative)
        logical,    intent(IN)  :: is_float     ! [--] Floating point or grounded?
        real(prec), intent(IN)  :: sigma(:)     ! nz   [--] Vertical sigma coordinates (sigma==height)
        real(prec), intent(IN)  :: sigt(:)      ! nz-1 [--] Vertical height axis (0:1) 
        real(prec), intent(IN)  :: dsigt_a(:)   ! nz-1 [--] d Vertical height axis (0:1) 
        real(prec), intent(IN)  :: dsigt_b(:)   ! nz-1 [--] d Vertical height axis (0:1) 
        real(prec), intent(IN)  :: dt           ! [a] Time step 

        ! Local variables 
        integer :: k, nz, nzt, nzz, ki  
        real(prec) :: T_pmp_base
        real(prec) :: Q_geo_now, ghf_conv  
        real(prec), allocatable :: cp_aa(:) 
        real(prec), allocatable :: ct_aa(:) 
        real(prec), allocatable :: advecxy_aa(:) 
        real(prec), allocatable :: advecz_aa(:) 
        real(prec), allocatable :: Q_strn_aa(:) 
        
        real(prec), allocatable :: subd(:)     ! nzz 
        real(prec), allocatable :: diag(:)     ! nzz  
        real(prec), allocatable :: supd(:)     ! nzz 
        real(prec), allocatable :: rhs(:)      ! nzz 
        real(prec), allocatable :: solution(:) ! nzz
        real(prec), allocatable :: factor(:)   ! nzt 
        real(prec) :: dsigmaBot
        
        nz  = size(sigma,1)
        nzt = size(sigt,1)   ! == nz-1, only layer midpoints where T is defined
        nzz = nz+1           ! == nz+1, layer midpoints plus upper and lower boundaries

        allocate(cp_aa(nzt))
        allocate(ct_aa(nzt))
        allocate(advecxy_aa(nzt))
        allocate(advecz_aa(nzt))
        allocate(Q_strn_aa(nzt))
        
        allocate(subd(nzz))
        allocate(diag(nzz))
        allocate(supd(nzz))
        allocate(rhs(nzz))
        allocate(solution(nzz))
        allocate(factor(nzt))

        ! Get average values over column for now 
        cp_aa = sum(cp) / size(cp,1)
        ct_aa = sum(ct) / size(ct,1)
        
        ! Get pressure melting point temperature at the base [celcius]
        T_pmp_base = calc_T_pmp(H_ice,sigma=0.0_prec,T0=0.0_prec)

        ! Get geothermal heat flux in proper units 
        Q_geo_now = Q_geo*1e-3*sec_year   ! [mW m-2] => [J m-2 a-1]

        ! Get advection on layer midpoints (aa nodes)
        do k = 1, nzt 
            advecxy_aa(k) = 0.5*(advecxy(k)+advecxy(k+1))
        end do 

        ! Get strain heating on layer midpoints (aa nodes)
        do k = 1, nzt 
            Q_strn_aa(k) = 0.5*(Q_strn(k)+Q_strn(k+1))
        end do 

        ! Step 1: apply vertical advection 
!         call advection_1D_upwind(advecz_aa,T_ice,uz,H_ice,T_srf,T_base,sigt)
        !call advection_1D_upwind_2ndorder(advecz_aa,T_ice,uz,H_ice,T_srf,T_base,sigt)
        call advection_1D_laxwendroff(advecz_aa,T_ice,uz,H_ice,T_srf,T_base,sigt,dt)
        !T_ice = T_ice - dt*advecz_aa 

        ! Step 2: apply vertical diffusion 

        ! Ice surface 
        subd(nzt+2) = 0.0_prec
        diag(nzt+2) = 1.0_prec
        supd(nzt+2) = 0.0_prec
        rhs(nzt+2)  = T_srf

        ! Ice interior, layers 1:nzt  (matrix elements 2:nzt+1)

        factor = dt * ct_aa / (rho_ice*cp_aa) / H_ice**2

        subd(2:nzt+1) = -factor * dsigt_a(1:nzt)
        supd(2:nzt+1) = -factor * dsigt_b(1:nzt)
        diag(2:nzt+1) = 1.0_prec - subd(2:nzt+1) - supd(2:nzt+1)
        rhs(2:nzt+1)  = T_ice(1:nzt) + Q_strn_aa(1:nzt)*dt - dt*advecz_aa 

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
                dsigmaBot = sigt(1) - 0.0 

                ! backward Euler flux basal boundary condition
                subd(1) =  0.0_prec
                diag(1) =  1.0_prec
                supd(1) = -1.0_prec
                rhs(1)  = (Q_b + Q_geo_now) * dsigmaBot*H_ice / ct_aa(1)

            end if   ! melting or frozen

        end if  ! floating or grounded 

        ! Call solver 
        call tridiag_solver(subd,diag,supd,solution,rhs)

        ! Copy the solution into the temperature variables
        !T_srf        = solution(1)
        T_ice(1:nzt) = solution(2:nzt+1)
        T_base       = solution(1)

        return 

    end subroutine mali_temp_diffusion_column

    subroutine advection_1D_upwind(advecz,Q,uz,H_ice,Q_srf,Q_base,sigt)
        ! Calculate vertical advection term advecz, which enters
        ! advection equation as
        ! Q_new = Q - dt*advecz = Q - dt*u*dQ/dx

        implicit none 

        real(prec), intent(OUT)   :: advecz(:) ! nzt, cell centers
        real(prec), intent(INOUT) :: Q(:)      ! nzt, cell centers
        real(prec), intent(IN)    :: uz(:)     ! nzt+1 == nz, cell boundaries
        real(prec), intent(IN)    :: H_ice     ! Ice thickness 
        real(prec), intent(IN)    :: Q_srf     ! Surface value
        real(prec), intent(IN)    :: Q_base    ! Base value
        real(prec), intent(IN)    :: sigt(:)   ! nzt, cell centers
        
        ! Local variables
        integer :: k, nzt   
        real(prec) :: u_aa 
        real(prec) :: dx 

        nzt = size(sigt,1)

        ! At center of base layer
        k    = 1 
        u_aa = 0.5*(uz(k)+uz(k+1)) ! Get velocity at cell center
        if (u_aa > 0.0) then 
            ! Upwind positive
            dx       = H_ice*(sigt(k+1)-sigt(k))
            advecz(k) = uz(k+1)*(Q(k+1)-Q(k))/dx   
        else
            ! Upwind negative
            dx       = H_ice*(sigt(k)-0.0)
            advecz(k) = uz(k)*(Q(k)-Q_base)/dx
        end if 

        ! Loop over internal cell centers and perform upwind advection 
        do k = 2, nzt-1 
            u_aa = 0.5*(uz(k)+uz(k+1)) ! Get velocity at cell center
            if (u_aa > 0.0) then 
                ! Upwind positive
                dx = H_ice*(sigt(k+1)-sigt(k))
                advecz(k) = uz(k+1)*(Q(k+1)-Q(k))/dx   
            else
                ! Upwind negative
                dx = H_ice*(sigt(k)-sigt(k-1))
                advecz(k) = uz(k)*(Q(k)-Q(k-1))/dx
            end if 
        end do 

        ! At center of surface layer
        k    = nzt 
        u_aa = 0.5*(uz(k)+uz(k+1)) ! Get velocity at cell center
        if (u_aa > 0.0) then 
            ! Upwind positive
            dx   = H_ice*(1.0-sigt(k))
            advecz(k) = uz(k+1)*(Q_srf-Q(k))/dx   
        else
            ! Upwind negative
            dx       = H_ice*(sigt(k)-sigt(k-1))
            advecz(k) = uz(k)*(Q(k)-Q(k-1))/dx
        end if 

        ! Replace old 
        return 

    end subroutine advection_1D_upwind

    subroutine advection_1D_upwind_2ndorder(advecz,Q,uz,H_ice,Q_srf,Q_base,sigt)
        ! Calculate vertical advection term advecz, which enters
        ! advection equation as
        ! Q_new = Q - dt*advecz = Q - dt*u*dQ/dx

        implicit none 

        real(prec), intent(OUT)   :: advecz(:) ! nzt, cell centers
        real(prec), intent(INOUT) :: Q(:)      ! nzt, cell centers
        real(prec), intent(IN)    :: uz(:)     ! nzt+1 == nz, cell boundaries
        real(prec), intent(IN)    :: H_ice     ! Ice thickness 
        real(prec), intent(IN)    :: Q_srf     ! Surface value
        real(prec), intent(IN)    :: Q_base    ! Base value
        real(prec), intent(IN)    :: sigt(:)   ! nzt, cell centers
        
        ! Local variables
        integer :: k, nzt   
        real(prec) :: u_aa 
        real(prec) :: dx 

        nzt = size(sigt,1)

        ! At center of base layer
        k    = 1 
        u_aa = 0.5*(uz(k)+uz(k+1)) ! Get velocity at cell center
        if (u_aa > 0.0) then 
            ! Upwind positive
            dx       = H_ice*(sigt(k+2)-sigt(k))
            advecz(k) = uz(k+1)*(-Q(k+2)+4.0*Q(k+1)-3.0*Q(k))/(dx)   
        else
            ! Upwind negative
            dx       = H_ice*(sigt(k)-0.0)
            advecz(k) = uz(k)*(Q(k)-Q_base)/dx
        end if 

        ! At center of second layer 
        k = 2 
        u_aa = 0.5*(uz(k)+uz(k+1)) ! Get velocity at cell center
        if (u_aa > 0.0) then 
            ! Upwind positive
            dx = H_ice*(sigt(k+2)-sigt(k))
            advecz(k) = uz(k+1)*(-Q(k+2)+4.0*Q(k+1)-3.0*Q(k))/(dx)    
        else
            ! Upwind negative
            dx = H_ice*(sigt(k)-sigt(k-1))
            advecz(k) = uz(k)*(Q(k)-Q(k-1))/dx
        end if 

        ! Loop over internal cell centers and perform upwind advection 
        do k = 3, nzt-2 
            u_aa = 0.5*(uz(k)+uz(k+1)) ! Get velocity at cell center
            if (u_aa > 0.0) then 
                ! Upwind positive
                dx = H_ice*(sigt(k+2)-sigt(k))
                advecz(k) = uz(k+1)*(-Q(k+2)+4.0*Q(k+1)-3.0*Q(k))/(dx)   
            else
                ! Upwind negative
                dx = H_ice*(sigt(k)-sigt(k-2))
                advecz(k) = uz(k)*(3.0*Q(k)-4.0*Q(k-1)+Q(k-2))/(dx)
            end if 
        end do 

        ! At center of layer below surface layer
        k = nzt-1  
        u_aa = 0.5*(uz(k)+uz(k+1)) ! Get velocity at cell center
        if (u_aa > 0.0) then 
            ! Upwind positive
            dx = H_ice*(sigt(k+1)-sigt(k))
            advecz(k) = uz(k+1)*(Q(k+1)-Q(k))/dx   
        else
            ! Upwind negative
            dx = H_ice*(sigt(k)-sigt(k-2))
            advecz(k) = uz(k)*(3.0*Q(k)-4.0*Q(k-1)+Q(k-2))/(dx)
        end if

        ! At center of surface layer
        k    = nzt 
        u_aa = 0.5*(uz(k)+uz(k+1)) ! Get velocity at cell center
        if (u_aa > 0.0) then 
            ! Upwind positive
            dx   = H_ice*(1.0-sigt(k))
            advecz(k) = uz(k+1)*(Q_srf-Q(k))/dx   
        else
            ! Upwind negative
            dx       = H_ice*(sigt(k)-sigt(k-2))
            advecz(k) = uz(k)*(3.0*Q(k)-4.0*Q(k-1)+Q(k-2))/(dx)
        end if 

        ! Replace old 
        return 

    end subroutine advection_1D_upwind_2ndorder

    subroutine advection_1D_laxwendroff(advecz,Q,uz,H_ice,Q_srf,Q_base,sigt,dt)
        ! Calculate vertical advection term advecz, which enters
        ! advection equation as
        ! Q_new = Q - dt*advecz = Q - dt*u*dQ/dx
        ! from: https://www.uni-muenster.de/imperia/md/content/physik_tp/lectures/ws2016-2017/num_methods_i/advection.pdf

        ! BROKEN 

        implicit none 

        real(prec), intent(OUT)   :: advecz(:) ! nzt, cell centers
        real(prec), intent(INOUT) :: Q(:)      ! nzt, cell centers
        real(prec), intent(IN)    :: uz(:)     ! nzt+1 == nz, cell boundaries
        real(prec), intent(IN)    :: H_ice     ! Ice thickness 
        real(prec), intent(IN)    :: Q_srf     ! Surface value
        real(prec), intent(IN)    :: Q_base    ! Base value
        real(prec), intent(IN)    :: sigt(:)   ! nzt, cell centers
        real(prec), intent(IN)    :: dt 

        ! Local variables
        integer :: k, nzt   
        real(prec) :: u_aa, u_acz, Q_lo, Q_hi  
        real(prec) :: dx 

        nzt = size(sigt,1)

        ! At center of base layer
        k    = 1 
        u_aa = 0.5*(uz(k)+uz(k+1)) ! Get velocity at cell center
        
!         dx   = sigt(k) - 0.0 
!         Q_lo = 0.5*(Q(k)+Q_base) - (uz(k)*dt/(2.0*dx))*(Q(k)-Q_base)
        Q_lo = Q_base 
        
        dx   = sigt(k) - 0.0 
        Q_hi = 0.5*(Q(k+1)+Q(k)) - (uz(k+1)*dt/(2.0*dx))*(Q(k+1)-Q(k))
        
        dx   = 0.5*(sigt(k)+sigt(k+1)) - 0.0        ! 0.5*(sigt(k)+sigt(k+1)) - 0.5*(sigt(k-1)+sigt(k))
        advecz(k) = u_aa*(Q_hi-Q_lo)/dx 

        ! Loop over internal cell centers
        do k = 2, nzt-1 
            
            u_aa = 0.5*(uz(k)+uz(k+1)) ! Get velocity at cell center
            
            dx   = sigt(k) - sigt(k-1) 
            Q_lo = 0.5*(Q(k)+Q(k-1)) - (uz(k)*dt/(2.0*dx))*(Q(k)-Q(k-1))

            dx   = sigt(k+1) - sigt(k)
            Q_hi = 0.5*(Q(k+1)+Q(k)) - (uz(k+1)*dt/(2.0*dx))*(Q(k+1)-Q(k))
            
            dx   = 0.5*(sigt(k+1) - sigt(k-1))     ! 0.5*(sigt(k)+sigt(k+1)) - 0.5*(sigt(k-1)+sigt(k)) 
            advecz(k) = u_aa*(Q_hi-Q_lo)/dx 

        end do 

        ! At center of surface layer
        k    = nzt 
        u_aa = 0.5*(uz(k)+uz(k+1)) ! Get velocity at cell center
        
        dx   = sigt(k) - sigt(k-1) 
        Q_lo = 0.5*(Q(k)+Q(k-1)) - (uz(k)*dt/(2.0*dx))*(Q(k)-Q(k-1))

!         dx   = 1.0 - sigt(k)
!         Q_hi = 0.5*(Q_srf+Q(k)) - (uz(k+1)*dt/(2.0*dx))*(Q_srf-Q(k))
        Q_hi = Q_srf 

        dx   = 1.0 - 0.5*(sigt(k-1)+sigt(k))        ! 0.5*(sigt(k)+sigt(k+1)) - 0.5*(sigt(k-1)+sigt(k))
        advecz(k) = u_aa*(Q_hi-Q_lo)/dx 
         
        return 

    end subroutine advection_1D_laxwendroff

        subroutine temperature_matrix_elements(&
         deltat,                &
         nVertLevels,           &
         layerCenterSigma,      &
         dsigmaTerm,            &
         iceCond_kt,            &
         floatingMask,          &
         thickness,             &
         temperature,           &
         surfaceTemperature,    &
         basalTemperature,      &
         heatDissipation,       &
         basalHeatFlux,         &
         basalFrictionFlux,     &
         subd,                  &
         diag,                  &
         supd,                  &
         rhs)

      ! Note: Matrix elements (subd, supd, diag, rhs) are indexed from 1 to nVertLevels+2,
      !       whereas temperature is indexed from 1 to nVertLevels.
      !       The first row of the matrix is the equation for surfaceTemperature, and
      !       the last row is the equation for basalTemperature.

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real(kind=prec), intent(in) :: &
           deltat                   !< Input: time step (s)

      integer, intent(in) :: &
           nVertLevels              !< Input: number of vertical layers

      real(kind=prec), dimension(nVertLevels), intent(in) :: &
           layerCentersigma         !< Input: sigma coordinate at layer midpoints

      real(kind=prec), dimension(:,:), intent(in) :: &
           dsigmaTerm               !< Input: vertical grid quantities

      real(kind=prec), dimension(nVertLevels), intent(in) :: &
           iceCond_kt               !< Input: vertical grid quantities

      integer, intent(in) :: &
           floatingMask             !< Input: = 1 where ice is floating, else = 0

      real(kind=prec), intent(in) ::  &
           thickness                !< Input: ice thickness (m)

      real(kind=prec), dimension(nVertLevels), intent(in) :: &
           temperature              !< Input: ice temperature (deg C)

      real(kind=prec), intent(in) :: &
           surfaceTemperature,    & !< Input: surface ice temperature (deg C)
           basalTemperature         !< Input: basal ice temperature (deg C)

      real(kind=prec), dimension(nVertLevels), intent(in) :: &
           heatDissipation          !< Input: interior heat dissipation (deg/s)

      real(kind=prec), intent(in) :: &
           basalHeatFlux            !< Input: geothermal heat flux (W m-2), positive up

      real(kind=prec), intent(in) :: &
           basalFrictionFlux        !< Input: basal friction heat flux (W m-2), >= 0

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real(kind=prec), dimension(:), intent(out) :: &
           subd, diag, supd, rhs    !< Output: tridiagonal matrix elements

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      real(kind=prec) :: pmpTemperatureBed  ! pressure melting temperature at bed
      real(kind=prec) :: factor             ! term in tridiagonal matrix
      real(kind=prec) :: dsigmaBot          ! bottom layer thickness in sigma coords

      real(kind=prec) :: iceConductivity   
      real(kind=prec) :: cp_ice   
      

      ! Get mean column ice conductivity 
      iceConductivity = sum(iceCond_kt) / size(iceCond_kt,1)


      ! Compute subdiagonal, diagonal, and superdiagonal matrix elements,
      ! along with the right-hand-side vector

      ! upper boundary: set to surfaceTemperature

      subd(1) = 0.0_prec
      diag(1) = 1.0_prec
      supd(1) = 0.0_prec
      rhs(1)  = surfaceTemperature

      ! ice interior, layers 1:nVertLevels  (matrix elements 2:nVertLevels+1)

      factor = deltat * iceConductivity / (rho_ice*cp_ice) / thickness**2

      subd(2:nVertLevels+1) = -factor * dsigmaTerm(1:nVertLevels,1)
      supd(2:nVertLevels+1) = -factor * dsigmaTerm(1:nVertLevels,2)
      diag(2:nVertLevels+1) = 1.0_prec - subd(2:nVertLevels+1) - supd(2:nVertLevels+1)
      rhs(2:nVertLevels+1)  = temperature(1:nVertLevels) + heatDissipation(1:nVertLevels)*deltat

      ! basal boundary:
      ! for grounded ice, a heat flux is applied
      ! for floating ice, the basal temperature is held constant

      !Note: If T(upn) < T_pmp, then require dT/dsigma = H/k * (G + taub*ubas)
      !       That is, net heat flux at lower boundary must equal zero.
      !      If T(upn) >= Tpmp, then set T(upn) = Tpmp

      if (floatingMask == 1) then

         subd(nVertLevels+2) = 0.0_prec
         diag(nVertLevels+2) = 1.0_prec
         supd(nVertLevels+2) = 0.0_prec
         rhs(nVertLevels+2)  = basalTemperature

      else    ! grounded ice

        ! ajr: commented out for compiling
         !call pressure_melting_point(thickness, pmpTemperatureBed)

         if (abs(basalTemperature - pmpTemperatureBed) < 0.001_prec) then

            ! hold basal temperature at pressure melting point

            subd(nVertLevels+2) = 0.0_prec
            diag(nVertLevels+2) = 1.0_prec
            supd(nVertLevels+2) = 0.0_prec
            rhs(nVertLevels+2)  = pmpTemperatureBed

         else   ! frozen at bed
                ! maintain balance of heat sources and sinks
                ! (conductive flux, geothermal flux, and basal friction)

            ! Note: basalHeatFlux is generally >= 0, since defined as positive up

            ! calculate dsigma for the bottom layer between the basal boundary and the temperature point above
            dsigmaBot = 1.0_prec - layerCentersigma(nVertLevels)

            ! backward Euler flux basal boundary condition
            subd(nVertLevels+2) = -1.0_prec
            diag(nVertLevels+2) =  1.0_prec
            supd(nVertLevels+2) =  0.0_prec
            rhs(nVertLevels+2)  = (basalFrictionFlux + basalHeatFlux) * dsigmaBot*thickness / iceConductivity

         endif   ! melting or frozen

      end if     ! floating or grounded

    end subroutine temperature_matrix_elements


    subroutine calc_sigt_terms(dsig_a,dsig_b,sigmamid,sigma)
        ! sigma = depth axis (1: base, nz: surface)
        ! Calculate ak, bk terms as defined in Hoffmann et al (2018)
        implicit none 

        real(prec), intent(INOUT) :: dsig_a(:)    ! nz-1
        real(prec), intent(INOUT) :: dsig_b(:)    ! nz-1
        real(prec), intent(INOUT) :: sigmamid(:)        ! nz-1 
        real(prec), intent(IN)    :: sigma(:)           ! nz 

        ! Local variables 
        integer :: k, nz_layers, nz  

        nz = size(sigma)
        nz_layers = nz - 1  

        ! Get sigmamid (midpoints of sigma layers - between sigma values)
        do k = 1, nz_layers
            sigmamid(k) = 0.5 * (sigma(k+1) + sigma(k))
        end do 


        k = 1 
        dsig_a(k) = 1.0 / ((sigma(k+1)-sigma(k))*(sigmamid(k)-sigma(k)))

        do k = 2, nz_layers
            dsig_a(k) = 1.0/ ( (sigma(k+1) - sigma(k)) * &
                                (sigmamid(k) - sigmamid(k-1)) )
        enddo

        do k = 1, nz_layers-1
            dsig_b(k) = 1.0/ ( (sigma(k+1) - sigma(k)) * &
                                (sigmamid(k+1) - sigmamid(k)) )
        end do

        k = nz_layers
        dsig_b(k) = 1.0/( (sigma(k+1) - sigma(k)) * &
                            (sigma(k+1) - sigmamid(k)) )

        return 

    end subroutine calc_sigt_terms

    subroutine tridiag_solver(a,b,c,x,y)
        !|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        !
        !  !  routine tridiag_solver
        !
        !> \brief MPAS solve tridiagonal matrix
        !> \author William Lipscomb
        !> \date   October 2015
        !> \details
        !>  This routine solves a tridiagonal matrix equation, given the matrix
        !>  coefficients and right-hand side.
        !-----------------------------------------------------------------------

        !-----------------------------------------------------------------
        ! input variables
        !-----------------------------------------------------------------

        real(kind=prec), dimension(:), intent(in)  :: a !< Input: Lower diagonal; a(1) is ignored
        real(kind=prec), dimension(:), intent(in)  :: b !< Input: Main diagonal
        real(kind=prec), dimension(:), intent(in)  :: c !< Input: Upper diagonal; c(n) is ignored
        real(kind=prec), dimension(:), intent(in)  :: y !< Input: Right-hand side

        !-----------------------------------------------------------------
        ! output variables
        !-----------------------------------------------------------------

        real(kind=prec), dimension(:), intent(out) :: x !< Output: Unknown vector

        !-----------------------------------------------------------------
        ! local variables
        !-----------------------------------------------------------------

        real(kind=prec), dimension(size(a)) :: aa
        real(kind=prec), dimension(size(a)) :: bb

        integer :: n, i

        n = size(a)

        aa(1) = c(1) / b(1)
        bb(1) = y(1) / b(1)

        do i = 2, n
            aa(i) = c(i) / (b(i)-a(i)*aa(i-1))
            bb(i) = (y(i)-a(i)*bb(i-1)) / (b(i)-a(i)*aa(i-1))
        end do

        x(n) = bb(n)

        do i = n-1, 1, -1
            x(i) = bb(i) - aa(i)*x(i+1)
        end do

        return 

    end subroutine tridiag_solver

end module icetemp_mali


