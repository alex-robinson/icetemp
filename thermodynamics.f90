module thermodynamics 
    ! This module contains some general thermodynamics subroutines
    ! that could be used by many icetemp solver approaches.
    ! Note: once icetemp is working well, this module could be 
    ! remerged into icetemp as one module. 

    use defs, only : prec, sec_year, pi, T0 

    implicit none 

    real(prec), parameter :: L_ice = 3.35e5       ! Specific latent heat of fusion of ice [J Kg-1]

    private 
    public :: calc_advec_horizontal_column
    public :: calc_basal_temp_gradients_column
    public :: calc_bmb_grounded
    public :: calc_strain_heating
    public :: calc_basal_heating
    public :: calc_specific_heat_capacity
    public :: calc_thermal_conductivity
    public :: calc_T_pmp
    public :: calc_f_pmp 
    public :: my_robin_solution 
    public :: error_function
    
contains 


        subroutine calc_advec_horizontal_column(advecxy,var_ice,ux,uy,dx,i,j)
        ! Newly implemented advection algorithms (ajr) 
        ! Output: [K a-1 m-2]

        implicit none

        real(prec), intent(OUT) :: advecxy(:) 
        real(prec), intent(IN)  :: var_ice(:,:,:)   ! Enth, T, age, etc...
        real(prec), intent(IN)  :: ux(:,:,:) 
        real(prec), intent(IN)  :: uy(:,:,:)
        real(prec), intent(IN)  :: dx  
        integer,    intent(IN)  :: i, j 

        ! Local variables 
        integer :: k, nx, ny, nz 
        real(prec) :: ux_aa, uy_aa 
        real(prec) :: dx_inv, dx_inv2
        real(prec) :: advecx, advecy 

        ! Define some constants 
        dx_inv  = 1.0_prec / dx 
        dx_inv2 = 1.0_prec / (2.0_prec*dx)

        nx = size(var_ice,1)
        ny = size(var_ice,2)
        nz = size(var_ice,3) 

        ! Loop over each point in the column
        do k = 1, nz 

            ! Estimate direction of current flow into cell (x and y)
            ux_aa = 0.5_prec*(ux(i,j,k)+ux(i-1,j,k))
            uy_aa = 0.5_prec*(uy(i,j,k)+uy(i,j-1,k))

            ! Explicit form (to test different order approximations)
            if (ux_aa .gt. 0.0 .and. i .ge. 3) then  
                ! Flow to the right 

                ! 1st order
!                 advecx = dx_inv * ux(i-1,j,k)*(var_ice(i,j,k)-var_ice(i-1,j,k))
                ! 2nd order
                advecx = dx_inv2 * ux(i-1,j,k)*(-(4.0*var_ice(i-1,j,k)-var_ice(i-2,j,k)-3.0*var_ice(i,j,k)))
            
            else if (ux_aa .lt. 0.0 .and. i .le. nx-2) then 
                ! Flow to the left

                ! 1st order 
!                 advecx = dx_inv * ux(i,j,k)*(var_ice(i+1,j,k)-var_ice(i,j,k))
                ! 2nd order
                advecx = dx_inv2 * ux(i,j,k)*((4.0*var_ice(i+1,j,k)-var_ice(i+2,j,k)-3.0*var_ice(i,j,k)))
            
            else 
                ! No flow 
                advecx = 0.0

            end if 

            if (uy_aa .gt. 0.0 .and. j .ge. 3) then   
                ! Flow to the right 

                ! 1st order
!                 advecy = dx_inv * uy(i,j-1,k)*(var_ice(i,j,k)-var_ice(i,j-1,k))
                ! 2nd order
                advecy = dx_inv2 * uy(i,j-1,k)*(-(4.0*var_ice(i,j-1,k)-var_ice(i,j-2,k)-3.0*var_ice(i,j,k)))
            
            else if (uy_aa .lt. 0.0 .and. j .le. ny-2) then 
                ! Flow to the left

                ! 1st order 
!                 advecy = dx_inv * uy(i,j,k)*(var_ice(i,j+1,k)-var_ice(i,j,k))
                ! 2nd order
                advecy = dx_inv2 * uy(i,j,k)*((4.0*var_ice(i,j+1,k)-var_ice(i,j+2,k)-3.0*var_ice(i,j,k)))
            
            else
                ! No flow 
                advecy = 0.0 

            end if 
                    
            ! Combine advection terms for total contribution 
            advecxy(k) = (advecx+advecy)

        end do 

        return 

    end subroutine calc_advec_horizontal_column

    subroutine calc_basal_temp_gradients_column(dTdz_b,dTrdz_b,T_ice,T_rock,H_ice,f_grnd,sigma,sigmar,H_rock)
        ! Calculate the temperature gradients in basal layer of ice and upper layer of bedrock 
        ! Returns dTdz_b and dTrdz_b in [K/m], positive upwards 

        real(prec), intent(OUT) :: dTdz_b 
        real(prec), intent(OUT) :: dTrdz_b 
        real(prec), intent(IN) :: T_ice(:) 
        real(prec), intent(IN) :: T_rock(:) 
        real(prec), intent(IN) :: H_ice 
        real(prec), intent(IN) :: f_grnd 
        real(prec), intent(IN) :: sigma(:) 
        real(prec), intent(IN) :: sigmar(:) 
        real(prec), intent(IN) :: H_rock 

        ! Local variables 
        integer    :: nzr  
        real(prec) :: dz 

        nzr = size(T_rock,1) 

        ! Get gradient in ice 
        if (H_ice .gt. 0.0_prec) then 
            dz = H_ice * (sigma(2)-sigma(1))
            dTdz_b = (T_ice(2) - T_ice(1)) / dz 
        else 
            dTdz_b = 0.0_prec 
        end if 

        ! Get gradient in rock 
        dz = H_rock * (sigmar(nzr)-sigma(nzr-1))
        dTrdz_b = (T_rock(nzr) - T_rock(nzr-1)) / dz  
        
        return 

    end subroutine calc_basal_temp_gradients_column

    elemental subroutine calc_bmb_grounded(bmb_grnd,T_prime_b,dTdz_b,dTrdz_b,kt_b,rho_ice, &
                                            Q_b,f_grnd,kt_m)
        ! Calculate everywhere there is at least some grounded ice 
        ! (centered aa node calculation)

        implicit none 
        
        real(prec), intent(OUT) :: bmb_grnd          ! [m/a] Basal mass balance, grounded
        real(prec), intent(IN)  :: T_prime_b         ! [K] Basal ice temp relative to pressure melting point (ie T_prime_b=0 K == temperate)
        real(prec), intent(IN)  :: dTdz_b            ! [K/m] Gradient of temperature in ice, basal layer
        real(prec), intent(IN)  :: dTrdz_b           ! [K/m] Gradient of temperature in rock, upper layer
        real(prec), intent(IN)  :: kt_b              ! [J a-1 m-1 K-1] Heat conductivity in ice, basal layer 
        real(prec), intent(IN)  :: rho_ice           ! [kg m-3] Ice density 
        real(prec), intent(IN)  :: Q_b               ! [J a-1 m-2] Basal heat production from friction and strain heating
        real(prec), intent(IN)  :: f_grnd            ! Grounded fraction (centered aa node)                 
        real(prec), intent(IN)  :: kt_m              ! [J a-1 m-1 K-1] Heat conductivity in mantle (lithosphere) 

        ! Local variables
        real(prec) :: Q_geo_now, coeff 
        real(prec), parameter :: tol = 1e-10 

        ! Get geothermal heat flux in the right units
!         Q_geo_now = Q_geo*1e-3*sec_year    ! [mW m-2] => [J a-1 m-2]

        ! Calculate the grounded basal mass balance following 
        ! Cuffey and Patterson (2010), Eq. 9.38 (Page 420)
        ! with an addition  of the term Q_strn_b 
        ! Note: formula only applies when base is temperate, ie
        ! when f_pmp > 0.0 
        if (T_prime_b .eq. 0.0_prec) then 
            ! Bed is temperate, calculate basal mass balance  

!                 if (cond_bed) then 
!                     ! Following grisli formulation: 
                    
!                     bmb_grnd = -1.0_prec/(rho_ice*L_ice)* ( Q_b + kt_b*dTdz_b - kt_m*dTrdz_b ) 

!                 else
!                     ! Classic Cuffey and Patterson (2010) formula 

!                     bmb_grnd = -1.0_prec/(rho_ice*L_ice)* ( Q_b + kt_b*dTdz_b + (Q_geo_now) ) 

!                 end if 
            
            bmb_grnd = -1.0_prec/(rho_ice*L_ice)* ( Q_b + kt_b*dTdz_b - kt_m*dTrdz_b ) 

        else 
            ! No basal mass change possible if bed is not temperate 

            bmb_grnd = 0.0_prec 

        end if 

        ! Limit small values to avoid underflow errors 
        if (abs(bmb_grnd) .lt. tol) bmb_grnd = 0.0_prec 

        return 

    end subroutine calc_bmb_grounded 

    subroutine calc_strain_heating(Q_strn,de,visc,cp,rho_ice)

        ! Calculate the general 3D internal strain heating
        ! as sum(D_ij*tau_ij)  (strain*stress)
        ! where stress has been calculated as stress_ij = 2*visc*strain
        ! Units: Q_strn = Q * 1/(cp*rho) = [J a-1 m-3] * [K m3 J-1] = [K a-1]

        implicit none

        real(prec), intent(OUT) :: Q_strn(:,:,:)          ! [Pa m a-1 ??] Heat production
        !type(stress_3D_class), intent(IN) :: strss        ! Stress tensor
        !type(strain_3D_class), intent(IN) :: strn         ! Strain rate tensor
        real(prec),            intent(IN) :: de(:,:,:)    ! [UNITS?] Effective strain rate 
        real(prec),            intent(IN) :: visc(:,:,:)  ! [UNITS?] Viscosity
        real(prec),            intent(IN) :: cp(:,:,:)    ! [J/kg/K] Specific heat capacity
        real(prec),            intent(IN) :: rho_ice      ! [kg m-3] Ice density 

        ! Local variables
        integer :: i, j, k, nx, ny, nz 
        real(prec), parameter :: Q_strn_max = 0.1         ! Check this limit!! 

        nx = size(Q_strn,1)
        ny = size(Q_strn,2)
        nz = size(Q_strn,3)

        ! Note: we use the simpler approach because in the shallow
        ! model, the stress rate is simply the strain rate squared

        ! Directly from Q_strn = tr(stress*strain)/cp
        ! (Cuffey and Patterson (2010) pag. 417, eq. 9.30)

!         Q_strn = ( strss%txx*strn%dxx &
!                  + strss%tyy*strn%dyy &
!                  + strss%tzz*strn%dzz &    ! this term is not available yet in the code, needs to be calculated
!              + 2.0*strss%txy*strn%dxy &
!              + 2.0*strss%txz*strn%dxz &
!              + 2.0*strss%tyz*strn%dyz ) * 1.0/(cp*rho_ice) 

        ! Simpler approach:
        ! Calculate Q_strn from effective strain rate and viscosity
        ! (Greve and Blatter (2009) eqs. 4.7 and 5.65): 
        !     Q_strn = tr(stress*strn)/cp = tr(2*visc*strn*strn)/cp = 2*visc*tr(strn*strn)/cp = 4*visc*de^2/cp
        !     with tr(strn*strn) = 2*de^2

        Q_strn = 4.0*visc*de**2 * 1.0/(cp*rho_ice) 

        ! Limit strain heating to reasonable values 
        where (Q_strn .gt. Q_strn_max) Q_strn = Q_strn_max

        return 

    end subroutine calc_strain_heating
    
    subroutine calc_basal_heating(Q_b,Q_strn_b,ux_base,uy_base,taub_acx,taub_acy,f_grnd_acx,f_grnd_acy)
        ! Qb [J a-1 m-2] == [m a-1] * [J m-3]

        implicit none 

        real(prec), intent(OUT) :: Q_b(:,:)               ! [J a-1 K-1] Basal heat production (friction + strain)
        real(prec), intent(IN)  :: Q_strn_b(:,:)          ! Basal heat production (strain heating only)
        real(prec), intent(IN)  :: ux_base(:,:)           ! Basal velocity, x-component (staggered x)
        real(prec), intent(IN)  :: uy_base(:,:)           ! Basal velocity, y-compenent (staggered y)
        real(prec), intent(IN)  :: taub_acx(:,:)          ! Basal friction (staggered x)
        real(prec), intent(IN)  :: taub_acy(:,:)          ! Basal friction (staggered y)
        real(prec), intent(IN)  :: f_grnd_acx(:,:)        ! Grounded fraction (staggered x)         
        real(prec), intent(IN)  :: f_grnd_acy(:,:)        ! Grounded fraction (staggered y)         
        
        ! Local variables
        integer    :: i, j, nx, ny 
        real(prec), allocatable :: Qb_acx(:,:), Qb_acy(:,:)
        
        nx = size(Q_b,1)
        ny = size(Q_b,2)

        ! Allocate staggered friction heat variables 
        allocate(Qb_acx(nx,ny))
        allocate(Qb_acy(nx,ny)) 

        ! Determine basal frictional heating values (staggered acx/acy nodes)
        ! (scaled by fraction of cell that is grounded)
        Qb_acx = f_grnd_acx * abs(ux_base*taub_acx)   ! [Pa m a-1] == [J a-1 m-2]
        Qb_acy = f_grnd_acy * abs(uy_base*taub_acy)   ! [Pa m a-1] == [J a-1 m-2]

        ! Initially set basal heating to zero everywhere 
        Q_b = 0.0  

        ! Get basal frictional heating on centered nodes (Aa grid)
        do j = 2, ny-1
        do i = 2, nx-1

            ! Average from Ac nodes to Aa node
            Q_b(i,j) = 0.25*(Qb_acx(i,j)+Qb_acx(i-1,j)+Qb_acy(i,j)+Qb_acy(i,j-1))

        end do
        end do

        ! Add strain heating in basal layer of ice sheet 
        Q_b = Q_b + Q_strn_b 

        return 

    end subroutine calc_basal_heating

    elemental function calc_specific_heat_capacity(T_ice) result(cp)

        implicit none 

        real(prec), intent(IN) :: T_ice  
        real(prec) :: cp 

        ! Specific heat capacity (Greve and Blatter, 2009, Eq. 4.39; Ritz, 1987)
        cp = (146.3 +7.253*T_ice)    ! [J kg-1 K-1]

        return 

    end function calc_specific_heat_capacity

    elemental function calc_thermal_conductivity(T_ice) result(ct)

        implicit none 

        real(prec), intent(IN) :: T_ice  
        real(prec) :: ct 

        ! Heat conductivity (Greve and Blatter, 2009, Eq. 4.37; Ritz, 1987)
        ct = 9.828*exp(-0.0057*T_ice)*sec_year  ! [W m-1 K-1 * sec_year] => [J m-1 K-1 a-1]

        return 

    end function calc_thermal_conductivity
    
    elemental function calc_T_pmp(H_ice,T0) result(T_pmp)
        ! Greve and Blatter (Chpt 4, pg 54), Eq. 4.13
        ! This gives the pressure-corrected melting point of ice
        ! where H_ice is the thickness of ice overlying the current point 
        
        implicit none 

        real(prec), intent(IN) :: H_ice  ! [m] Ice thickness above this point
        real(prec), intent(IN) :: T0     ! [K] Reference freezing point of water (273.15 K)
        real(prec) :: T_pmp              ! [K] Pressure corrected melting point

        ! Local variables  
        real(prec), parameter :: beta1 = 8.74e-4    ! [K m^-1]   beta1 = (beta*rho*g), beta=9.8e-8 [K Pa^-1]

        T_pmp = T0 - beta1*H_ice 

        ! ajr: note: should we account here for whether ice is floating or not, changing the pressure? 
        
        return 

    end function calc_T_pmp 

    elemental function calc_f_pmp(T_ice,T_pmp,gamma,is_float) result(f_pmp)
        ! Calculate the fraction of gridpoint at the pressure melting point (pmp),
        ! ie, when T_ice >= T_pmp. Facilitates a smooth transition between
        ! frozen and temperate ice. (Greve, 2005; Hindmarsh and Le Meur, 2001)

        implicit none 

        real(prec), intent(IN) :: T_ice
        real(prec), intent(IN) :: T_pmp
        real(prec), intent(IN) :: gamma
        logical, intent(IN) :: is_float  
        real(prec) :: f_pmp 

        if (is_float) then
            ! Floating points are temperate by default
            f_pmp = 1.0 

        else 
            ! Calculate the fraction at the pressure melting point 

            if (gamma .eq. 0.0) then
                ! No decay function, binary pmp fraction

                if (T_ice .ge. T_pmp) then 
                    f_pmp = 1.0
                else 
                    f_pmp = 0.0 
                end if 

            else

                ! Apply decay function 
                f_pmp = min(1.0, exp((T_ice-T_pmp)/gamma) )

            end if 

        end if 

        ! Ensure pure zero values below a threshold 
        f_pmp = max(f_pmp, 1e-4)

        return 

    end function calc_f_pmp 
    
    function my_robin_solution(sigma,T_pmp,kt,cp,rho_ice,H_ice,T_srf,mb_net,Q_geo,is_float) result(T_ice)
        ! This function will impose a temperature solution in a given ice column.
        ! For:
        !  Grounded ice with positive net mass balance: Robin solution where possible
        !  Grounded ice with negative net mass balance: Linear profile 
        !  Floating ice: Linear profile 
        !  No or thin ice: Surface temperature 
        
        implicit none 

        real(prec), intent(IN) :: sigma(:) 
        real(prec), intent(IN) :: T_pmp(:)
        real(prec), intent(IN) :: kt(:) 
        real(prec), intent(IN) :: cp(:) 
        real(prec), intent(IN) :: rho_ice 
        real(prec), intent(IN) :: H_ice 
        real(prec), intent(IN) :: T_srf 
        real(prec), intent(IN) :: mb_net
        real(prec), intent(IN) :: Q_geo 
        logical,    intent(IN) :: is_float 

        real(prec) :: T_ice(size(sigma,1))

        ! Local variables 
        integer    :: k, nz 
        real(prec) :: dTdz_b, z, kappa, ll   

        real(prec), parameter :: sqrt_pi   = sqrt(pi) 
        real(prec), parameter :: T_ocn     = 271.15   ! [K]
        real(prec), parameter :: H_ice_min = 0.1      ! [m] Minimum ice thickness to calculate Robin solution 
        real(prec), parameter :: mb_net_min = 1e-2    ! [m a-1] Minimum allowed net mass balance for stability
        real(prec) :: Q_geo_now, mb_now  

        nz = size(T_ice,1) 
        
        Q_geo_now = Q_geo *1e-3*sec_year    ! [mW m-2] => [J a-1 m-2]

        ! Calculate temperature gradient at base 
        dTdz_b = -Q_geo_now/kt(1) 

        if (.not. is_float .and. H_ice .gt. H_ice_min .and. mb_net .gt. 0.0) then 
            ! Impose Robin solution 
            
            !mb_now = max(mb_net,mb_net_min)
            mb_now = mb_net 

            do k = 1, nz 
                z     = sigma(k)*H_ice              ! Ice thickness up to this layer from base 
                kappa = kt(k)/(cp(k)*rho_ice)       ! Thermal diffusivity 
                ll    = sqrt(2*kappa*H_ice/mb_now)  ! Thermal_length_scale

                ! Calculate ice temperature for this layer 
                T_ice(k) = (sqrt_pi/2.0)*ll*dTdz_b*(error_function(z/ll)-error_function(H_ice/ll)) + T_srf 
            end do 

        else if (.not. is_float .and. H_ice .gt. H_ice_min) then 
            ! Impose linear profile with temperate base

            T_ice(nz) = T_srf
            T_ice(1)  = T_pmp(1)

            ! Intermediate layers are linearly interpolated 
            do k = 2, nz-1 
                T_ice(k) = T_ice(1)+sigma(k)*(T_ice(nz)-T_ice(1))
            end do 
            
        else if (is_float) then 
            ! Floating - impose linear profile between sea and surface temps 

            T_ice(nz) = T_srf
            T_ice(1)  = T_ocn 

            ! Intermediate layers are linearly interpolated 
            do k = 2, nz-1 
                T_ice(k) = T_ice(1)+sigma(k)*(T_ice(nz)-T_ice(1))
            end do 
            
        else 
            ! Ice thickness too small 
            T_ice = T_srf 

        end if 

        return 

    end function my_robin_solution

    function error_function(X) result(ERR)
        ! Purpose: Compute error function erf(x)
        ! Input:   x   --- Argument of erf(x)
        ! Output:  ERR --- erf(x)
        
        implicit none 

        real(prec), intent(IN)  :: X
        real(prec) :: ERR
        
        ! Local variables:
        real(prec)              :: EPS
        real(prec)              :: X2
        real(prec)              :: ER
        real(prec)              :: R
        real(prec)              :: C0
        integer                 :: k
        
        EPS = 1.0e-15
        X2  = X * X
        if (abs(X) < 3.5) then
            ER = 1.0
            R  = 1.0
            do k = 1, 50
                R  = R * X2 / (real(k, prec) + 0.5)
                ER = ER+R
                if(abs(R) < abs(ER) * EPS) then
                    C0  = 2.0 / sqrt(pi) * X * exp(-X2)
                    ERR = C0 * ER
                    EXIT
                end if
            end do
        else
            ER = 1.0
            R  = 1.0
            do k = 1, 12
                R  = -R * (real(k, prec) - 0.5) / X2
                ER = ER + R
                C0  = EXP(-X2) / (abs(X) * sqrt(pi))
                ERR = 1.0 - C0 * ER
                if(X < 0.0) ERR = -ERR
            end do
        end if

        return

    end function error_function

end module thermodynamics 
