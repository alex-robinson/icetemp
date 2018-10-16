module icetemp_grisli
    ! Wrapping the original grisli icetemp solution for yelmo 

    use defs, only : prec, pi, g, sec_year, T0, rho_ice, rho_sw, rho_w 
    use solver_tridiagonal, only : tridiag 
    use thermodynamics, only : calc_advec_horizontal_column_sico1, calc_advec_horizontal_column

    implicit none
    
    ! Impose parameter choice here for testing conductive bedrock or not
    logical,    parameter :: conductive_bedrock = .FALSE. 
    
    private
    public :: calc_icetemp_grisli_3D_up
    public :: calc_icetemp_grisli_column_up

contains 
    
    subroutine calc_icetemp_grisli_3D_up(T_ice,T_rock,T_pmp,cp,ct,ux,uy,uz,Q_strn,Q_b, &
                                            Q_geo,T_srf,H_ice,H_w,smb,bmb,f_grnd,sigma,dt,dx)
        ! GRISLI solver for thermodynamics for a given column of ice 
        ! Note sigma=height, k=1 base, k=nz surface 
        ! Note T_ice, T_pmp in [degC], not [K]

        implicit none 

        real(prec), intent(OUT) :: T_ice(:,:,:)     ! [degC] Ice column temperature
        real(prec), intent(OUT) :: T_rock(:,:,:)    ! [degC] Bedrock column temperature
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
        real(prec), intent(IN)  :: sigma(:)         ! [--] Vertical sigma coordinates (sigma==height)
        real(prec), intent(IN)  :: dt               ! [a] Time step 
        real(prec), intent(IN)  :: dx               ! [a] Horizontal grid step 

        ! Local variable
        integer :: i, j, k, nx, ny, nz  
        real(prec), allocatable  :: advecxy(:)   ! [K a-1 m-2] Horizontal heat advection 
        logical :: is_float 
        real(prec) :: H_w_dot 

        nx = size(T_ice,1)
        ny = size(T_ice,2)
        nz = size(T_ice,3)

        allocate(advecxy(nz))
        advecxy = 0.0 

        do j = 3, ny-2
        do i = 3, nx-2 

            ! Calculate the contribution of horizontal advection to column solution
            call calc_advec_horizontal_column(advecxy,T_ice,ux,uy,dx,i,j)
!             call calc_advec_horizontal_column_sico1(advecxy,T_ice,ux,uy,dx,i,j)

!             if (H_ice(16,16) .gt. 2800.0 .and. i .eq. 24 .and. j .eq. 16) then
!                 do k = 1, nz  
!                     write(*,*) advecxy(k), Q_strn(i,j,k)
!                 end do 
!                 stop 
!             end if 
            
            is_float = (f_grnd(i,j) .eq. 0.0)
             
            ! Call thermodynamic solver for the column (solver expects degrees Celcius)

            T_ice(i,j,:)  = T_ice(i,j,:)  - T0 
            T_rock(i,j,:) = T_rock(i,j,:) - T0

            call calc_icetemp_grisli_column_up(T_ice(i,j,:),T_rock(i,j,:),T_pmp(i,j,:)-T0,cp(i,j,:),ct(i,j,:), &
                                            uz(i,j,:),Q_strn(i,j,:),advecxy,Q_b(i,j),Q_geo(i,j),T_srf(i,j)-T0, &
                                            H_ice(i,j),H_w(i,j),bmb(i,j),is_float,sigma,dt)
            
            T_ice(i,j,:)  = T_ice(i,j,:)  + T0 
            T_rock(i,j,:) = T_rock(i,j,:) + T0 
            
        end do 
        end do 

        return 

    end subroutine calc_icetemp_grisli_3D_up

    subroutine calc_icetemp_grisli_column_up(T_ice,T_rock,T_pmp,cp,ct,uz,Q_strn,advecxy,Q_b, &
                                            Q_geo,T_srf,H_ice,H_w,bmb,is_float,sigma,dt)
        ! GRISLI solver for thermodynamics for a given column of ice 
        ! Note sigma=height, k=1 base, k=nz surface 
        ! Note T_ice, T_pmp in [degC], not [K]

        implicit none 

        real(prec), intent(OUT) :: T_ice(:)     ! [degC] Ice column temperature
        real(prec), intent(OUT) :: T_rock(:)    ! [degC] Bedrock column temperature
        real(prec), intent(IN)  :: T_pmp(:)     ! [degC] Pressure melting point temp.
        real(prec), intent(IN)  :: cp(:)        ! [J kg-1 K-1] Specific heat capacity
        real(prec), intent(IN)  :: ct(:)        ! [J a-1 m-1 K-1] Heat conductivity 
        real(prec), intent(IN)  :: uz(:)        ! [m a-1] Vertical velocity 
        real(prec), intent(IN)  :: Q_strn(:)    ! [K a-1] Internal strain heat production in ice
        real(prec), intent(IN)  :: advecxy(:)   ! [K a-1] Horizontal heat advection 
        real(prec), intent(IN)  :: Q_b          ! [J a-1 m-2] Basal frictional heat production 
        real(prec), intent(IN)  :: Q_geo        ! [mW m-2] Geothermal heat flux 
        real(prec), intent(IN)  :: T_srf        ! [degC] Surface temperature 
        real(prec), intent(IN)  :: H_ice        ! [m] Ice thickness 
        real(prec), intent(IN)  :: H_w          ! [m] Basal water layer thickness 
        real(prec), intent(IN)  :: bmb          ! [m a-1] Basal mass balance (melting is negative)
        logical,    intent(IN)  :: is_float     ! [--] Floating point or grounded?
        real(prec), intent(IN)  :: sigma(:)     ! [--] Vertical sigma coordinates (sigma==height)
        real(prec), intent(IN)  :: dt           ! [a] Time step 

        ! Local variables 
        integer    :: k, ki, nz   
        integer    :: nzm, nzz, nzb   
        real(prec) :: dzm, de, da, ro, rom, romg, row, cl, cm, cpm
        real(prec) :: dzz, dah, dou, dzi  
        real(prec) :: ct_bas, ct_haut 
        real(prec) :: dah_bas, dzz_bas, dah_haut, dzz_haut
        real(prec) :: ctm 
        real(prec) :: acof1, bcof1, ccof1, s0mer, tbmer 
        real(prec) :: tbdot, tdot, tss, tbb
        real(prec) :: bmelt
        real(prec) :: Q_geo_now
        real(prec) :: H_w_dot    
        integer    :: ifail 

        real(prec), allocatable :: aa(:), bb(:), cc(:), rr(:), hh(:) 
        real(prec), allocatable :: T(:), T_new(:) 

        real(prec), allocatable :: advecz(:) 

        ! Store local vector sizes 
        nz      = size(T_ice,1)    ! Number of ice points 
        nzm     = size(T_rock,1)   ! Number of bedrock points 
        nzz     = nz+nzm           ! Total grid points (ice plus bedrock)

        nzb     = nzm + 1          ! Index of ice base in column (ice plus bedrock)

        ! Check that sigma is equally spaced
        de = sigma(2)-sigma(1)
        do k = 2, nz 
            dzi = sigma(k)-sigma(k-1)
            if (abs(de-dzi).gt.1e-2) then 
                write(*,*) "calc_icetemp:: Error: currently only evenly spaced vertical axis is allowed."
                write(*,*) "sigma = ", sigma 
                stop 
            end if 
        end do 

        ! Some parameters defined here for now, for testing 
        ! (these will move to a parameter file later)
        dzm     = 600.0            ! [m] Bedrock step height 
        da      = 4.0e7            ! Mantle diffusion
        ro      = rho_ice          ! [kg m-3] Ice density 
        rom     = 3300.0           ! [kg m-3] Density of mantle
        romg    = rom*g            ! [kg m-2 s-2] Density * gravity constant
        row     = rho_sw           ! [kg m-3] Seawater density 
        cl      = 3.35e5           ! [J kg-1] specific latent heat of fusion of ice
        cm      = 1.04e8           ! [J m-1 K-1 a-1] Thermal conductivity of mantle 
        cpm     = 1000.0           ! [J kg-1 K-1] Specific heat capacity of mantle

        ! Coefficients for sea temperature from Jenkins (1991)
        acof1   = -0.0575
        bcof1   =  0.0901
        ccof1   = 7.61e-4
        s0mer   = 34.75

        ! Convert units of Geothermal heat flux
        Q_geo_now = Q_geo *sec_year*1e-3   ! [mW/m2] => [J a-1 m-2]

        ! Derived constants
        ctm     = dt*cm/rom/cpm/dzm/dzm  

        ! Allocate local variables 
        allocate(aa(nzz),bb(nzz),cc(nzz),rr(nzz),hh(nzz))
        allocate(T(nzz),T_new(nzz))

        allocate(advecz(nz))


        !call advection_1D_upwind_interp(advecz,T_ice(2:nz-1),uz,H_ice,T_ice(nz),T_ice(1),sigma(2:nz-1),sigma)

        ! Populate local temperature column (ice+rock) 
        T(1:nzm)       = T_rock
        T(nzb:nzz)     = T_ice 

        ! Populate new solution to be safe 
        T_new = T 

        ! Calculate the sea temperature below ice shelf if point is floating 
        if (is_float) then
            tbmer = acof1*s0mer + bcof1 + ccof1*H_ice*ro/row
        else 
            tbmer = 0.0 
        end if 

        ! Determine expected change in basal water total for this timestep [m]
        H_w_dot = -(bmb*(rho_w/rho_ice))*dt 

        ! Bedrock conditions ===========

        if (conductive_bedrock) then 
            ! With conductive bedrock 

            ! Boundary conditions at the base of the bedrock 
            aa(1) = -1.0
            bb(1) =  1.0
            cc(1) =  0.0
            rr(1) = dzm*Q_geo_now/cm
            
            ! Internal bedrock points 
            do k = 2, nzm
                aa(k) = -ctm
                bb(k) = 1.0+2.0*ctm
                cc(k) = -ctm
                rr(k) = T(k)
            end do
            
        else 
            ! Without conductive bedrock: fill in bedrock solution 

            ! Prescribe solution for bedrock points:
            ! Calculate temperature in the bedrock as a linear gradient of ghf
            ! Note: Using T(nzb) from the previous timestep for simplicity
            do k = 1, nzm
                aa(k) = 0.0 
                bb(k) = 1.0 
                cc(k) = 0.0 
                rr(k) = T(nzb)+dzm*(nzm-k+1)*Q_geo_now/cm
            end do
            
        end if 

        ! Ice sheet conditions ===========

        if (H_ice .gt. 10.0) then
            ! General case (H_ice>10m)
            
            ! Limits at the base of the ice sheet 

            ! Frozen base 
            if ( (.not. is_float) .and.  &
                 (T(nzb).lt.(T_pmp(1)) .or. H_w+H_w_dot .lt. 0.0_prec) ) then 
                ! Frozen, grounded bed; or about to become so with removal of basal water

                de = (sigma(2)-sigma(1))*H_ice 

                if (conductive_bedrock) then
                    ! With conductive bedrock 
                    dzi     = de*cm/ct(1)
                    aa(nzb) = -dzm/(dzm+dzi)
                    bb(nzb) = 1.0
                    cc(nzb) = -dzi/(dzm+dzi)
                    rr(nzb) = de*Q_b/ct(1)*dzm/(dzm+dzi)

                else
                    ! Without conductive bedrock 
                    aa(nzb) =  0.0
                    bb(nzb) =  1.0
                    cc(nzb) = -1.0
                    rr(nzb) = (Q_geo_now+Q_b)/ct(1)*de

                end if

            else
                ! Temperate ice or shelf ice 

                aa(nzb) = 0.0
                bb(nzb) = 1.0
                cc(nzb) = 0.0

                if (.not. is_float) then
                    rr(nzb) = T_pmp(1)
                else
                    rr(nzb) = Tbmer
                end if

            end if

            ! Internal ice sheet points
            
            ct_haut = 2.0*(ct(1)*ct(2))/(ct(1)+ct(2))
!             ct_haut = calc_harmonic_mean_wtd(ct(1),ct(2),sigma(2)-sigma(1),sigma(3)-sigma(2))

            do k = nzb+1, nzz-1

                ki = k - nzm     ! Ice sheet index 

!                 de      = H_ice*(sigma(ki)-sigma(ki-1))
!                 dah_bas = dt/de
!                 dzz_bas = dt/(de*de) * 1.0/(cp(ki)*ro)

!                 de       = H_ice*(sigma(ki+1)-sigma(ki))
!                 dah_haut = dt/de
!                 dzz_haut = dt/(de*de) * 1.0/(cp(ki)*ro)

!                 dzz = ((sigma(ki)-sigma(ki-1))*dzz_bas+(sigma(ki+1)-sigma(ki))*dzz_haut) &
!                         / (sigma(ki+1)-sigma(ki-1))

                de  = H_ice*(sigma(ki+1)-sigma(ki)) 
                dah = dt/de
                dzz = dt/(de*de) * 1.0/(cp(ki)*ro)

                dah_bas  = dah 
                dah_haut = dah 
                dzz_bas  = dzz 
                dzz_haut = dzz 

                ! Conductivity as the Harmonic average of the grid points
                ct_bas  = ct_haut
                ct_haut = 2.0*(ct(ki)*ct(ki+1))/(ct(ki)+ct(ki+1))
!                 ct_haut = calc_harmonic_mean_wtd(ct(ki),ct(ki+1),sigma(ki)-sigma(ki-1),sigma(ki+1)-sigma(ki))
            
                ! Vertical advection (centered)
                aa(k) = -dzz_bas*ct_bas   - uz(ki)*dah_bas/2.0
                bb(k) = 1.0+dzz*(ct_bas+ct_haut)
                cc(k) = -dzz_haut*ct_haut + uz(ki)*dah_haut/2.0
                rr(k) = T(k) + dt*Q_strn(ki) - dt*advecxy(ki) !- dt*advecz(ki)

            end do
            
            ! Prescribe surface temperature as boundary condition

            T(nzz)  = min(0.0,T_srf)

            aa(nzz) = 0.0
            bb(nzz) = 1.0
            cc(nzz) = 0.0
            rr(nzz) = T(nzz)
            
        else 
            ! Thin or no ice (H_ice<=10m)

            if (H_ice .gt. 0.0) then
                ! Ice exists, impose the gradient from the bedrock 

                ! Impose surface temperature at the surface
                tss = min(0.0,T_srf)

                if (.not. is_float .and. conductive_bedrock) then 
                    ! Grounded point with conductive bedrock 

                    ! Slope is based on gradient between bedrock and ice point
!                     de  = (sigma(2)-sigma(1))*H_ice 
!                     dou = ( (T(nzm)-T(nzb))/dzm*cm) / ct(1)*de
!                     tbb = tss+dou*(nz-1)
                    tbb = T(nzb) 

                else if (.not. is_float) then 
                    ! Grounded point 

                    ! Slope is based on geothermal heat flux directly 
                    tbb = tss + Q_geo_now/ct(1)*H_ice

                else 
                    ! Floating point 

                    ! Basal temperature is ocean temperature
                    tbb = tbmer
                     
                end if 

                ! Interpolate intermediate points
                do k = nzb, nzz
                    ki    = k-nzm  
                    aa(k) = 0.0
                    bb(k) = 1.0 
                    cc(k) = 0.0 
                    rr(k) = tbb+(tss-tbb)*sigma(ki)
                end do
                    
            else
                ! Ice-free point - set all points equal to pressure-melting point
                ! (ie, T0)

                aa(nzb:nzz) = 0.0
                bb(nzb:nzz) = 1.0 
                cc(nzb:nzz) = 0.0
                rr(nzb:nzz) = T_pmp

            end if

        end if 

        ! Solve temperature equation 
        ! =============================================================

        call tridiag(aa,bb,cc,rr,hh,nzz,ifail)

        ! Check tridiag results
        if (ifail .eq. 1) then
            write(*,*) 'Error in tridiag solver.'
            write(*,*) 'a, b, c, r, u'
            do k = 1, nzz
                write(*,*)'k = ', k
                write(*,*) aa(k),bb(k),cc(k),rr(k),hh(k)
            end do
            stop
        end if

        ! Fill in new temperature solution 
        T_new = hh 
        
        ! Populate T vector
        ! Apply limits: less than 3degC / 5 yrs variation and no colder than -70 degC
        do k = 1, nzz

            if (H_ice .gt. 10.0) then

                tdot = (T_new(k)-T(k))/dt

                if (k .lt. nzz .and. tdot .lt. -3.0) tdot = -3.0
                if (k .lt. nzz .and. tdot .gt.  3.0) tdot =  3.0

                T_new(k) = T(k) + tdot

                if (T_new(k) .lt. -70.0) T_new(k) = -70.0

            end if

        end do

        ! Diagnose rate of temperature change at ice base 
        tbdot = (T_new(nzb)-T(nzb))/dt

        ! Limit ice temperatures to the pressure melting point
        if (H_ice .gt. 0.0) then  
            do k = nzb, nzz
                ki = k - nzm 
                if (T_new(k) .gt. T_pmp(ki)) T_new(k) = T_pmp(ki)
            end do
        end if 

        ! Store solution back in output variables 
        T_rock = T_new(1:nzm)   
        T_ice  = T_new(nzb:nzz)

        return 

    end subroutine calc_icetemp_grisli_column_up 

    function calc_harmonic_mean(k1,k2) result(kmean)
        ! From wikipedia:
        ! https://en.wikipedia.org/wiki/Harmonic_mean

        implicit none 

        real(prec), intent(IN) :: k1, k2
        real(prec) :: kmean 

        kmean = 2.0*(k1*k2) / (k1+k2)

        return 

    end function calc_harmonic_mean
    
    function calc_harmonic_mean_wtd(k1,k2,w1,w2) result(kmean)
        ! From wikipedia:
        ! https://en.wikipedia.org/wiki/Harmonic_mean#Weighted_harmonic_mean

        implicit none 

        real(prec), intent(IN) :: k1, k2, w1, w2 
        real(prec) :: kmean 

        kmean = (w1+w2) / (w1/k1 + w2/k2)

        return 

    end function calc_harmonic_mean_wtd


    subroutine calc_icetemp_grisli_column_up_evensteps(T_ice,T_rock,T_pmp,cp,ct,uz,Q_strn,advecxy,Q_b, &
                                            Q_geo,T_srf,H_ice,H_w,bmb,is_float,sigma,dt)
        ! GRISLI solver for thermodynamics for a given column of ice 
        ! Note sigma=height, k=1 base, k=nz surface 
        ! Note T_ice, T_pmp in [degC], not [K]

        implicit none 

        real(prec), intent(OUT) :: T_ice(:)     ! [degC] Ice column temperature
        real(prec), intent(OUT) :: T_rock(:)    ! [degC] Bedrock column temperature
        real(prec), intent(IN)  :: T_pmp(:)     ! [degC] Pressure melting point temp.
        real(prec), intent(IN)  :: cp(:)        ! [J kg-1 K-1] Specific heat capacity
        real(prec), intent(IN)  :: ct(:)        ! [J a-1 m-1 K-1] Heat conductivity 
        real(prec), intent(IN)  :: uz(:)        ! [m a-1] Vertical velocity 
        real(prec), intent(IN)  :: Q_strn(:)    ! [K a-1] Internal strain heat production in ice
        real(prec), intent(IN)  :: advecxy(:)   ! [K a-1 m-2] Horizontal heat advection 
        real(prec), intent(IN)  :: Q_b          ! [J a-1 m-2] Basal frictional heat production 
        real(prec), intent(IN)  :: Q_geo        ! [mW m-2] Geothermal heat flux 
        real(prec), intent(IN)  :: T_srf        ! [degC] Surface temperature 
        real(prec), intent(IN)  :: H_ice        ! [m] Ice thickness 
        real(prec), intent(IN)  :: H_w          ! [m] Basal water layer thickness 
        real(prec), intent(IN)  :: bmb          ! [m a-1] Basal mass balance (melting is negative)
        logical,    intent(IN)  :: is_float     ! [--] Floating point or grounded?
        real(prec), intent(IN)  :: sigma(:)     ! [--] Vertical sigma coordinates (sigma==height)
        real(prec), intent(IN)  :: dt           ! [a] Time step 

        ! Local variables 
        integer    :: k, ki, nz   
        integer    :: nzm, nzz, nzb   
        real(prec) :: dzm, de, da, ro, rom, romg, row, cl, cm, cpm
        real(prec) :: dzz, dah, dou, dzi  
        real(prec) :: ct_bas, ct_haut 
        real(prec) :: ctm 
        real(prec) :: acof1, bcof1, ccof1, s0mer, tbmer 
        real(prec) :: tbdot, tdot, tss
        real(prec) :: bmelt
        real(prec) :: Q_geo_now
        real(prec) :: H_w_dot    
        integer    :: ifail 

        real(prec), allocatable :: aa(:), bb(:), cc(:), rr(:), hh(:) 
        real(prec), allocatable :: abis(:), bbis(:), cbis(:), rbis(:), hbis(:)  
        real(prec), allocatable :: T(:), T_new(:) 

        ! Store local vector sizes 
        nz      = size(T_ice,1)    ! Number of ice points 
        nzm     = size(T_rock,1)   ! Number of bedrock points 
        nzz     = nz+nzm           ! Total grid points (ice plus bedrock)

        nzb     = nzm + 1          ! Index of ice base in column (ice plus bedrock)


        ! Some parameters defined here for now, for testing 
        ! (these will move to a parameter file later)
        dzm     = 600.0            ! [m] Bedrock step height 
        de      = 1.0/(nz-1)       ! vertical step in ice
        da      = 4.0e7            ! Mantle diffusion
        ro      = rho_ice          ! [kg m-3] Ice density 
        rom     = 3300.0           ! [kg m-3] Density of mantle
        romg    = rom*g            ! [kg m-2 s-2] Density * gravity constant
        row     = rho_sw           ! [kg m-3] Seawater density 
        cl      = 3.35e5           ! [J kg-1] specific latent heat of fusion of ice
        cm      = 1.04e8           ! [J m-1 K-1 a-1] Thermal conductivity of mantle 
        cpm     = 1000.0           ! [J kg-1 K-1] Specific heat capacity of mantle

        ! Coefficients for sea temperature from Jenkins (1991)
        acof1   = -0.0575
        bcof1   =  0.0901
        ccof1   = 7.61e-4
        s0mer   = 34.75

        ! Convert units of Geothermal heat flux
        Q_geo_now = Q_geo *sec_year*1e-3   ! [mW/m2] => [J a-1 m-2]

        ! Derived constants
        ctm     = dt*cm/rom/cpm/dzm/dzm  

        ! Allocate local variables 
        allocate(aa(nzz),bb(nzz),cc(nzz),rr(nzz),hh(nzz))
        allocate(abis(nzz),bbis(nzz),cbis(nzz),rbis(nzz),hbis(nzz))

        allocate(T(nzz),T_new(nzz))

        ! Populate local temperature column (ice+rock) 
        T(1:nzm)       = T_rock 
        T(nzb:nzz)     = T_ice 

        ! Populate new solution to be safe 
        T_new = T 

        ! NOTE: Thermal properties are now calculated outside of the routine
        ! The below notes are included for reference. Also now cp is expected
        ! in units of [J kg-1 K-1] to be consistent with literature, while
        ! grisli originally used units cp*ro with units of [J m-3 K-1]

        ! GRISLI THERMAL PROPERTIES ===========================================
!         ! Calculate thermal properties
!         do k = 1, nz  
!             cp(k)    = (2115.3+7.79293*T_ice(k))             ! [J kg-1 K-1]
!             ct(k)    = 3.1014e8*exp(-0.0057*(T_ice(k)+T0))   ! [J m-1 K-1 a-1] 
!         end do 
        ! =====================================================================

        ! EISMINT THERMAL PROPERTIES ==========================================
!         cp = 2009.0        ! [J kg-1 K-1]
!         ct = 6.6e7         ! [J m-1 K-1 a-1] 
        ! =====================================================================
        
        ! Calculate the sea temperature below ice shelf if point is floating 
        if (is_float) then
            tbmer = acof1*s0mer + bcof1 + ccof1*H_ice*ro/row
        else 
            tbmer = 0.0 
        end if 

!         ! Diagnose basal melt rate and basal state (temperate, frozen, etc.)
!         call bmelt_grounded_column_up(bmelt,ibase,T_ice,T_rock,ct,Q_b,H_ice,H_w, &
!                                       Q_geo_now,is_float,de,dzm,cm,ro,cl,dt)


        ! Determine generated/lost water total for this timestep [m]
        H_w_dot = -(bmb*(rho_w/rho_ice))*dt 

        ! Bedrock conditions ===========

        if (conductive_bedrock) then 
            ! With conductive bedrock 

            ! Boundary conditions at the base of the bedrock 
            aa(1) = -1.0
            bb(1) =  1.0
            cc(1) =  0.0
            rr(1) = dzm*Q_geo_now/cm
            
            ! Internal bedrock points 
            do k = 2, nzm
                aa(k) = -ctm
                bb(k) = 1.0+2.0*ctm
                cc(k) = -ctm
                rr(k) = T(k)
            end do
            
        else 
            ! Without conductive bedrock: fill in bedrock solution 

            ! Prescribe solution for bedrock points:
            ! Calculate temperature in the bedrock as a linear gradient of ghf
            ! Note: Using T(nzb) from the previous timestep for simplicity
            do k = 1, nzm
                aa(k) = 0.0 
                bb(k) = 1.0 
                cc(k) = 0.0 
                rr(k) = T(nzb)+dzm*(nzm-k+1)*Q_geo_now/cm
            end do
            
        end if 

        ! Ice sheet conditions ===========

        if (H_ice .gt. 10.0) then
            ! General case (H_ice>10m)
            
            ! Limits at the base of the ice sheet 

            ! Frozen base 
!             if (.not. is_float .and. ((ibase.eq.1).or.(ibase.eq.4) &
!                 .or. ((ibase.eq.5).and.(T(nzb).lt.T_pmp(1)))) ) then
            
            if ( (.not. is_float) .and.  &
                 (T(nzb).lt.T_pmp(1) .or. H_w+H_w_dot .lt. 0.0_prec) ) then 
                ! Frozen, grounded bed; or about to become so with removal of basal water

                if (conductive_bedrock) then
                    ! With conductive bedrock 
                    dzi    = H_ice*de*cm/ct(1)
                    aa(nzb) = -dzm/(dzm+dzi)
                    bb(nzb) = 1.0
                    cc(nzb) = -dzi/(dzm+dzi)
                    rr(nzb) = H_ice*de*Q_b/ct(1)*dzm/(dzm+dzi)

                else
                    ! Without conductive bedrock 
                    aa(nzb) =  0.0
                    bb(nzb) =  1.0
                    cc(nzb) = -1.0
                    rr(nzb) = (Q_geo_now+Q_b)/ct(1)*H_ice*de

                end if

            else
                ! Temperate ice or shelf ice 

                aa(nzb) = 0.0
                bb(nzb) = 1.0
                cc(nzb) = 0.0

                if (.not. is_float) then
                    rr(nzb) = T_pmp(1)
                else
                    rr(nzb) = Tbmer
                end if

            end if

            ! Internal ice sheet points
        
            dou    = dt/de/de/H_ice/H_ice
            dah    = dt/H_ice/de
            ct_haut = 2.0*(ct(1)*ct(2))/(ct(1)+ct(2))

            do k = nzb+1, nzz-1

                ki = k - nzm     ! Ice sheet index 

                dzz = dou/(cp(ki)*ro)

                ! Conductivity as the Harmonic average of the grid points
                ct_bas  = ct_haut
                ct_haut = 2.0*(ct(ki)*ct(ki+1))/(ct(ki)+ct(ki+1))

                ! Vertical advection (centered)
                aa(k) = -dzz*ct_bas-uz(ki)*dah/2.0
                bb(k) = 1.0+dzz*(ct_bas+ct_haut)
                cc(k) = -dzz*ct_haut+uz(ki)*dah/2.0
                rr(k) = T(k)+dt*Q_strn(ki)-dt*advecxy(ki)

            end do
            
            ! Prescribe surface temperature as boundary condition

            T(nzz)  = min(0.0,T_srf)

            aa(nzz) = 0.0
            bb(nzz) = 1.0
            cc(nzz) = 0.0
            rr(nzz) = T(nzz)
            
        else 
            ! Thin or no ice (H_ice<=10m)

            if (H_ice .gt. 0.0) then
                ! Ice exists, impose the gradient from the bedrock 

                ! Impose surface temperature at the surface
                tss = min(0.0,T_srf)

                if (.not. is_float .and. conductive_bedrock) then 
                    ! Grounded point with conductive bedrock 

                    ! Slope is based on gradient between bedrock and ice point
                    dou = ( (T(nzm)-T(nzb))/dzm*cm) / ct(1)*de*H_ice

                else if (.not. is_float) then 
                    ! Grounded point 

                    ! Slope is based on geothermal heat flux directly 
                    dou = Q_geo_now/ct(1)*de*H_ice
                
                else 
                    ! Floating point 

                    ! Slope is gradient between surface and ocean temperature
                    dou = (tbmer-tss)*de
                     
                end if 

                ! Interpolate intermediate points
                do k = nzb, nzz 
                    aa(k) = 0.0
                    bb(k) = 1.0 
                    cc(k) = 0.0 
                    rr(k) = tss+dou*(nzz-k)
                end do
                    
            else
                ! Ice-free point

                aa(nzb:nzz) = 0.0
                bb(nzb:nzz) = 1.0 
                cc(nzb:nzz) = 0.0

                if (is_float) then 
                    rr(nzb:nzz) = tbmer
                else 
                    rr(nzb:nzz) = T_srf
                end if 

            end if

        end if 

        ! Solve temperature equation 
        ! =============================================================

        call tridiag(aa,bb,cc,rr,hh,nzz,ifail)

        ! Check tridiag results
        if (ifail .eq. 1) then
            write(*,*) 'Error in tridiag solver.'
            write(*,*) 'a, b, c, r, u'
            do k = 1, nzz
                write(*,*)'k = ', k
                write(*,*) aa(k),bb(k),cc(k),rr(k),hh(k)
            end do
            stop
        end if

        ! Fill in new temperature solution 
        T_new = hh 
        
        ! Diagnose rate of temperature change at ice base 
        tbdot = (T_new(1)-T(1))/dt

        ! Populate T vector
        ! Apply limits: less than 3degC / 5 yrs variation and no colder than -70 degC
        do k = 1, nzz

            if (H_ice .gt. 10.0) then

                tdot = (T_new(k)-T(k))/dt

                if (k .lt. nzz .and. tdot .lt. -3.0) tdot = -3.0
                if (k .lt. nzz .and. tdot .gt.  3.0) tdot =  3.0

                T_new(k) = T(k) + tdot

                if (T_new(k) .lt. -70.0) T_new(k) = -70.0

            end if

        end do

        ! Limit ice temperatures to the pressure melting point 
        do k = nzb, nzz
            ki = k - nzm 
            if (T_new(k) .gt. T_pmp(ki)) T_new(k) = T_pmp(ki)
        end do

        ! Store solution back in output variables 
        T_rock = T_new(1:nzm)
        T_ice  = T_new(nzb:nzz)

        return 

    end subroutine calc_icetemp_grisli_column_up_evensteps 

    subroutine advection_1D_upwind_interp(advecz,Q,uz,H_ice,Q_srf,Q_base,sigt,sigma)
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
        real(prec), intent(IN)    :: sigma(:)  ! nz (==nzt+1), cell borders
        
        ! Local variables
        integer :: k, nzt   
        real(prec) :: u_aa 
        real(prec) :: dx 
        
        real(prec), allocatable :: x1(:), Q1(:), uz1(:), advecz1(:)  
        integer :: nz1 
        real(prec) :: dx1 


        nzt = size(sigt,1)

        nz1 = 21 
        allocate(x1(nz1))
        allocate(Q1(nz1))
        allocate(uz1(nz1)) 
        allocate(advecz1(nz1)) 

        ! Generate hi-res axis 
        do k = 1, nz1 
            x1(k) = 0.0 + real(k-1,prec)* 1.0/real(nz1-1,prec)
        end do 

        dx1 = H_ice* (x1(2) - x1(1))

        ! Get interpolated values 
!         Q1  = interp_linear(x=[0.0_prec,sigt,1.0_prec],y=[Q_base,Q,Q_srf],xout=x1)
!         uz1 = interp_linear(x=sigma,y=uz,xout=x1)
        Q1  = interp_spline(x=[0.0_prec,sigt,1.0_prec],y=[Q_base,Q,Q_srf],xout=x1)
        uz1 = interp_spline(x=sigma,y=uz,xout=x1)

        ! At bottom of base layer
        k = 1 
        advecz1(k) = 0.0 

        ! Loop over internal cell centers and perform upwind advection 
        do k = 2, nz1-1 
            u_aa = uz1(k)  ! Get velocity at cell center
            if (u_aa > 0.0) then 
                ! Upwind positive
                advecz1(k) = 0.5*(uz1(k)+uz1(k+1))*(Q1(k+1)-Q1(k))/dx1   
            else
                ! Upwind negative
                advecz1(k) = 0.5*(uz1(k-1)+uz1(k))*(Q1(k)-Q1(k-1))/dx1
            end if 
        end do 

        ! At top of surface layer
        k = nz1 
        advecz1(k) = 0.0 

        ! Reinterpolate to old axis 
!         advecz = interp_linear(x=x1,y=advecz1,xout=sigt)
        advecz = interp_spline(x=x1,y=advecz1,xout=sigma)
        
        return 

    end subroutine advection_1D_upwind_interp
    
    ! 1D interpolation
    function interp_spline(x,y,xout) result(yout)

        implicit none 
 
        real(prec), dimension(:), intent(IN) :: x, y
        real(prec), dimension(:), intent(IN) :: xout
        real(prec), dimension(size(xout)) :: yout 
        real(prec), dimension(:), allocatable :: b, c, d 
        real(prec) :: uh, dx, yh  
        integer :: i, n, nout 

        n    = size(x) 
        nout = size(xout)

        ! Get spline coefficients b, c, d
        allocate(b(n),c(n),d(n))
        call spline (x, y, b, c, d, n)

        do i = 1, nout 
            if (xout(i) .lt. x(1)) then
                dx = x(1)-xout(i)
                uh = x(1)+dx
                yh = ispline(uh,x,y,b,c,d,n)
                yout(i) = y(1) + (y(1)-yh)
                !write(*,*) x(1), xout(i), dx, uh, y(1), yh, yout(i)
            else if (xout(i) .gt. x(n)) then
                dx = xout(i)-x(n)
                uh = x(n)-dx
                yh = ispline(uh,x,y,b,c,d,n)
                yout(i) = y(n) + (y(n)-yh)
                !write(*,*) x(n), xout(i), dx, uh, y(n), yh, yout(i)
            else
                yout(i) = ispline(xout(i), x, y, b, c, d, n)
            end if 
        end do 

        return

    contains 

        subroutine spline (x, y, b, c, d, n)
            !======================================================================
            !  SOURCE: http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch01/spline.f90
            !  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
            !  for cubic spline interpolation
            !  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
            !  for  x(i) <= x <= x(i+1)
            !  Alex G: January 2010
            !----------------------------------------------------------------------
            !  input..
            !  x = the arrays of data abscissas (in strictly increasing order)
            !  y = the arrays of data ordinates
            !  n = size of the arrays xi() and yi() (n>=2)
            !  output..
            !  b, c, d  = arrays of spline coefficients
            !  comments ...
            !  spline.f90 program is based on fortran version of program spline.f
            !  the accompanying function fspline can be used for interpolation
            !======================================================================
            implicit none
            integer n
            real(prec), dimension(:) :: x, y, b, c, d 
            !real(prec) x(n), y(n), b(n), c(n), d(n)
            integer i, j, gap
            real(prec) h

            gap = n-1
            ! check input
            if ( n < 2 ) return
            if ( n < 3 ) then
              b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
              c(1) = 0.
              d(1) = 0.
              b(2) = b(1)
              c(2) = 0.
              d(2) = 0.
              return
            end if
            !
            ! step 1: preparation
            !
            d(1) = x(2) - x(1)
            c(2) = (y(2) - y(1))/d(1)
            do i = 2, gap
              d(i) = x(i+1) - x(i)
              b(i) = 2.0*(d(i-1) + d(i))
              c(i+1) = (y(i+1) - y(i))/d(i)
              c(i) = c(i+1) - c(i)
            end do
            !
            ! step 2: end conditions 
            !
            b(1) = -d(1)
            b(n) = -d(n-1)
            c(1) = 0.0
            c(n) = 0.0
            if(n /= 3) then
              c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
              c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
              c(1) = c(1)*d(1)**2/(x(4)-x(1))
              c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
            end if
            !
            ! step 3: forward elimination 
            !
            do i = 2, n
              h = d(i-1)/b(i-1)
              b(i) = b(i) - h*d(i-1)
              c(i) = c(i) - h*c(i-1)
            end do
            !
            ! step 4: back substitution
            !
            c(n) = c(n)/b(n)
            do j = 1, gap
              i = n-j
              c(i) = (c(i) - d(i)*c(i+1))/b(i)
            end do
            !
            ! step 5: compute spline coefficients
            !
            b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
            do i = 1, gap
              b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
              d(i) = (c(i+1) - c(i))/d(i)
              c(i) = 3.*c(i)
            end do
            c(n) = 3.0*c(n)
            d(n) = d(n-1)
        end subroutine spline

        function ispline(u, x, y, b, c, d, n)
            !======================================================================
            ! SOURCE: http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch01/spline.f90
            ! function ispline evaluates the cubic spline interpolation at point z
            ! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
            ! where  x(i) <= u <= x(i+1)
            !----------------------------------------------------------------------
            ! input..
            ! u       = the abscissa at which the spline is to be evaluated
            ! x, y    = the arrays of given data points
            ! b, c, d = arrays of spline coefficients computed by spline
            ! n       = the number of data points
            ! output:
            ! ispline = interpolated value at point u
            !=======================================================================
            implicit none
            real(prec) ispline
            integer n
            ! real(prec)  u, x(n), y(n), b(n), c(n), d(n)
            real(prec) :: u 
            real(prec), dimension(:) :: x, y, b, c, d 
            integer i, j, k
            real(prec) dx

            ! if u is ouside the x() interval take a boundary value (left or right)
            if(u <= x(1)) then
              ispline = y(1)
              return
            end if
            if(u >= x(n)) then
              ispline = y(n)
              return
            end if

            !*
            !  binary search for for i, such that x(i) <= u <= x(i+1)
            !*
            i = 1
            j = n+1
            do while (j > i+1)
              k = (i+j)/2
              if(u < x(k)) then
                j=k
                else
                i=k
               end if
            end do
            !*
            !  evaluate spline interpolation
            !*
            dx = u - x(i)
            ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
        end function ispline

    end function interp_spline 

    function interp_linear(x,y,xout) result(yout)
        ! Interpolate y from ordered x to ordered xout positions

        implicit none 
 
        real(prec), dimension(:), intent(IN) :: x, y
        real(prec), dimension(:), intent(IN) :: xout
        real(prec), dimension(size(xout)) :: yout 
        integer :: i, j, n, nout 

        n    = size(x) 
        nout = size(xout)

!         write(*,*) minval(x), maxval(x), n, nout

        do i = 1, nout 
            if (xout(i) .lt. x(1)) then
                yout(i) = y(1)
!                 write(*,*) 1, xout(i)
            else if (xout(i) .gt. x(n)) then
                yout(i) = y(n)
!                 write(*,*) 2, xout(i)
            else
                do j = 1, n 
                    if (x(j) .ge. xout(i)) exit 
                end do

                if (j .eq. 1) then 
                    yout(i) = y(1) 
!                     write(*,*) 3, xout(i)
                else if (j .eq. n+1) then 
                    yout(i) = y(n)
!                     write(*,*) 4, xout(i)
                else 
                    yout(i) = interp_linear_internal(x(j-1:j),y(j-1:j),xout(i))
!                     write(*,*) 5, xout(i)
                end if 
            end if 
        end do

        return 

    contains 

            ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !   Subroutine :  interp_linear_internal
        !   Author     :  Alex Robinson
        !   Purpose    :  Interpolates for the y value at the desired x value, 
        !                 given x and y values around the desired point.
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function interp_linear_internal(x,y,xout) result(yout)

            implicit none

            real(prec), intent(IN)  :: x(2), y(2), xout
            real(prec) :: yout
            real(prec) :: alph

            if ( xout .lt. x(1) .or. xout .gt. x(2) ) then
                write(*,*) "interp1: xout < x0 or xout > x1 !"
                write(*,*) "xout = ",xout
                write(*,*) "x0   = ",x(1)
                write(*,*) "x1   = ",x(2)
                stop
            end if

            alph = (xout - x(1)) / (x(2) - x(1))
            yout = y(1) + alph*(y(2) - y(1))

            return

        end function interp_linear_internal 

    end function interp_linear

end module icetemp_grisli 

