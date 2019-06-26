module icetemp 
    ! Module contains the ice temperature and basal mass balance (grounded) solution

    use defs, only : prec, pi, g, sec_year, T0, rho_ice, rho_sw, rho_w, L_ice  
    use solver_tridiagonal, only : solve_tridiag 
    use thermodynamics, only : calc_bmb_grounded, calc_advec_vertical_column

    implicit none
    
    private
    public :: calc_temp_column 
    public :: calc_dzeta_terms

contains 

    subroutine calc_temp_column(T_ice,bmb_grnd,Q_ice_b,T_pmp,cp,kt,uz,Q_strn,advecxy,Q_b,Q_geo, &
                    T_srf,T_shlf,H_ice,H_w,f_grnd,zeta_aa,zeta_ac,dzeta_a,dzeta_b,dt)
        ! Thermodynamics solver for a given column of ice 
        ! Note zeta=height, k=1 base, k=nz surface 
        ! Note: nz = number of vertical boundaries (including zeta=0.0 and zeta=1.0), 
        ! temperature is defined for cell centers, plus a value at the surface and the base
        ! so nz_ac = nz_aa - 1 

        ! For notes on implicit form of advection terms, see eg http://farside.ph.utexas.edu/teaching/329/lectures/node90.html
        
        implicit none 

        real(prec), intent(INOUT) :: T_ice(:)     ! nz_aa [K] Ice column temperature
        real(prec), intent(INOUT) :: bmb_grnd     ! [m a-1] Basal mass balance (melting is negative)
        real(prec), intent(OUT)   :: Q_ice_b      ! [J a-1 m-2] Basal heat flux into ice (positive up)
        real(prec), intent(IN)    :: T_pmp(:)     ! nz_aa [K] Pressure melting point temp.
        real(prec), intent(IN)    :: cp(:)        ! nz_aa [J kg-1 K-1] Specific heat capacity
        real(prec), intent(IN)    :: kt(:)        ! nz_aa [J a-1 m-1 K-1] Heat conductivity 
        real(prec), intent(IN)    :: uz(:)        ! nz_ac [m a-1] Vertical velocity 
        real(prec), intent(IN)    :: Q_strn(:)    ! nz_aa [J a-1 m-3] Internal strain heat production in ice
        real(prec), intent(IN)    :: advecxy(:)   ! nz_aa [K a-1] Horizontal heat advection 
        real(prec), intent(IN)    :: Q_b          ! [J a-1 m-2] Basal frictional heat production 
        real(prec), intent(IN)    :: Q_geo        ! [mW m-2] Geothermal heat flux 
        real(prec), intent(IN)    :: T_srf        ! [K] Surface temperature 
        real(prec), intent(IN)    :: T_shlf       ! [K] Marine-shelf interface temperature
        real(prec), intent(IN)    :: H_ice        ! [m] Ice thickness 
        real(prec), intent(IN)    :: H_w          ! [m] Basal water layer thickness 
        real(prec), intent(IN)    :: f_grnd       ! [--] Grounded fraction
        real(prec), intent(IN)    :: zeta_aa(:)  ! nz_aa [--] Vertical sigma coordinates (zeta==height), layer centered aa-nodes
        real(prec), intent(IN)    :: zeta_ac(:)  ! nz_ac [--] Vertical height axis temperature (0:1), layer edges ac-nodes
        real(prec), intent(IN)    :: dzeta_a(:)  ! nz_aa [--] Solver discretization helper variable ak
        real(prec), intent(IN)    :: dzeta_b(:)  ! nz_aa [--] Solver discretization helper variable bk
        real(prec), intent(IN)    :: dt           ! [a] Time step 

        ! Local variables 
        integer :: k, nz_aa, nz_ac
        real(prec) :: Q_geo_now, ghf_conv 
        real(prec) :: Q_strn_now
        real(prec) :: H_w_predicted
        real(prec) :: T_excess
        real(prec) :: melt_internal   

        real(prec), allocatable :: advecz(:)   ! nz_aa, for explicit vertical advection solving
        logical, parameter      :: test_expl_advecz = .FALSE. 

        real(prec), allocatable :: kappa_aa(:)
        real(prec), allocatable :: dkappadz(:)

        real(prec), allocatable :: subd(:)     ! nz_aa 
        real(prec), allocatable :: diag(:)     ! nz_aa  
        real(prec), allocatable :: supd(:)     ! nz_aa 
        real(prec), allocatable :: rhs(:)      ! nz_aa 
        real(prec), allocatable :: solution(:) ! nz_aa
        real(prec) :: fac, fac_a, fac_b, uz_aa, dzeta, dz, dz1, dz2  
        real(prec) :: kappa_a, kappa_b, dza, dzb 

        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1)

        allocate(kappa_aa(nz_aa))
        allocate(dkappadz(nz_aa))
        
        allocate(subd(nz_aa))
        allocate(diag(nz_aa))
        allocate(supd(nz_aa))
        allocate(rhs(nz_aa))
        allocate(solution(nz_aa))

        ! Get geothermal heat flux in proper units 
        Q_geo_now = Q_geo*1e-3*sec_year   ! [mW m-2] => [J m-2 a-1]

        ! Calculate diffusivity on cell centers (aa-nodes)
        kappa_aa = kt / (rho_ice*cp)

        ! Calculate gradient in kappa (centered on aa-nodes)
        dkappadz = 0.0 
        do k = 2, nz_aa-1 
            dkappadz(k) = (kappa_aa(k+1)-kappa_aa(k-1))/(zeta_aa(k+1)-zeta_aa(k-1))
        end do 

        ! Step 1: apply vertical advection (for explicit testing)
        if (test_expl_advecz) then 
            allocate(advecz(nz_aa))
            advecz = 0.0
            call calc_advec_vertical_column(advecz,T_ice,uz,H_ice,zeta_aa)
            T_ice = T_ice - dt*advecz 
        end if 

        ! Step 2: apply vertical implicit diffusion-advection 
        
        ! == Ice base ==

        if (f_grnd .lt. 1.0) then
            ! Floating or partially floating ice - set temperature equal 
            ! to basal temperature at pressure melting point, or marine freezing temp,
            ! or weighted average between the two.

            ! Impose the weighted average of the pressure melting point and the marine freezing temp.
            subd(1) = 0.0_prec
            diag(1) = 1.0_prec
            supd(1) = 0.0_prec
            rhs(1)  = f_grnd*T_pmp(1) + (1.0-f_grnd)*T_shlf

        else 
            ! Grounded ice 

            ! Determine expected basal water thickness [m] for this timestep,
            ! using basal mass balance from previous time step (rough guess)
            H_w_predicted = H_w - (bmb_grnd*(rho_w/rho_ice))*dt 

            ! == Assign grounded basal boundary conditions ==

            if (T_ice(1) .lt. T_pmp(1) .or. H_w_predicted .lt. 0.0_prec) then   
                ! Frozen at bed, or about to become frozen 

                ! maintain balance of heat sources and sinks
                ! (conductive flux, geothermal flux, and basal friction)

                ! Note: basalHeatFlux is generally >= 0, since defined as positive up

                ! calculate dzeta for the bottom layer between the basal boundary and the temperature point above
                dzeta = zeta_aa(2) - zeta_aa(1)

                ! backward Euler flux basal boundary condition
                subd(1) =  0.0_prec
                diag(1) =  1.0_prec
                supd(1) = -1.0_prec
                rhs(1)  = (Q_b + Q_geo_now) * dzeta*H_ice / kt(1)
                
                if (.FALSE.) then 
                    ! Alternative boundary condition approach - works, but needs testing and refinement 

                ! Convert units of Q_strn [J a-1 m-3] => [K a-1]
                Q_strn_now = Q_strn(1)/(rho_ice*cp(1))

                fac      = dt * kt(1) / (rho_ice*cp(1)) / H_ice**2
               
                subd(1) =  0.0_prec
                diag(1) =  1.0_prec  + fac/((zeta_ac(2)-zeta_ac(1))*(zeta_aa(2) - zeta_aa(1))  )     !mmr  1.0_prec
                supd(1) =  -1.0_prec  * fac/((zeta_ac(2)-zeta_ac(1))*(zeta_aa(2) - zeta_aa(1))  )     !mmr -1.0_prec
                rhs(1)  =  T_ice(1) + fac* ((Q_geo_now) * H_ice / kt(1)) * 1.0_prec/( (zeta_ac(2)-zeta_ac(1))) + Q_strn_now*dt   !+ uz(1) * (Q_b_now) *H*dt / (( zeta_aa(2)-zeta_aa(1) ) * kt(1) )    ! mmr (Q_b_now) * dzetaBot*H / kt(1)
                
                end if 

            else 
                ! Temperate at bed 
                ! Hold basal temperature at pressure melting point

                subd(1) = 0.0_prec
                diag(1) = 1.0_prec
                supd(1) = 0.0_prec
                rhs(1)  = T_pmp(1) 

            end if   ! melting or frozen

        end if  ! floating or grounded 

        ! == Ice interior layers 2:nz_aa-1 ==

        do k = 2, nz_aa-1

            if (test_expl_advecz) then 
                ! No implicit vertical advection (diffusion only)
                uz_aa = 0.0 

            else
                ! With implicit vertical advection (diffusion + advection)
                uz_aa   = 0.5*(uz(k-1)+uz(k))   ! ac => aa nodes

            end if 

            ! Add 'vertical advection' due to the gradient in kappa 
            !uz_aa = uz_aa + dkappadz(k)

            ! Convert units of Q_strn [J a-1 m-3] => [K a-1]
            Q_strn_now = Q_strn(k)/(rho_ice*cp(k))

        if (.FALSE.) then 
            ! Stagger kappa to the lower and upper ac-nodes

            ! ac-node between k-1 and k 
            if (k .eq. 2) then 
                ! Bottom layer, kappa is kappa for now (later with bedrock kappa?)
                kappa_a = kappa_aa(1)
            else 
                ! Weighted average between lower half and upper half of point k-1 to k 
                dz1 = zeta_ac(k-1)-zeta_aa(k-1)
                dz2 = zeta_aa(k)-zeta_ac(k-1)
                kappa_a = (dz1*kappa_aa(k-1) + dz2*kappa_aa(k))/(dz1+dz2)
            end if 

            ! ac-node between k and k+1 

            ! Weighted average between lower half and upper half of point k to k+1
            dz1 = zeta_ac(k+1)-zeta_aa(k)
            dz2 = zeta_aa(k+1)-zeta_ac(k+1)
            kappa_b = (dz1*kappa_aa(k) + dz2*kappa_aa(k+1))/(dz1+dz2)

        else 
            ! ajr: simply use aa-node kappas for now
            kappa_a = kappa_aa(k) 
            kappa_b = kappa_a

        end if 

            ! Vertical distance for centered difference advection scheme
            dz      =  H_ice*(zeta_aa(k+1)-zeta_aa(k-1))
            
            fac_a   = -kappa_a*dzeta_a(k)*dt/H_ice**2
            fac_b   = -kappa_b*dzeta_b(k)*dt/H_ice**2

            subd(k) = fac_a - uz_aa * dt/dz
            supd(k) = fac_b + uz_aa * dt/dz
            diag(k) = 1.0_prec - fac_a - fac_b
            rhs(k)  = T_ice(k) + dt*Q_strn_now - dt*advecxy(k) 

        end do 

        ! == Ice surface ==

        subd(nz_aa) = 0.0_prec
        diag(nz_aa) = 1.0_prec
        supd(nz_aa) = 0.0_prec
        rhs(nz_aa)  = min(T_srf,T0)


        ! == Call solver ==

        call solve_tridiag(subd,diag,supd,rhs,solution)

        ! Copy the solution into the temperature variables
        T_ice  = solution
        
        ! === Treat basal mass balance and high temperatures ===
        
        ! Calculate heat flux at ice base as temperature gradient * conductivity [J a-1 m-2]
        if (H_ice .gt. 0.0_prec) then 
            dz = H_ice * (zeta_aa(2)-zeta_aa(1))
            Q_ice_b = kt(1) * (T_ice(2) - T_ice(1)) / dz 
        else 
            Q_ice_b = 0.0  
        end if 
        
        ! First calculate internal melt (only allow melting, no accretion)
        
        melt_internal = 0.0 

        do k = nz_aa-1, 2, -1 
            ! Descend from surface to base layer (center of layer)

            ! Store temperature difference above pressure melting point (excess energy)
            T_excess = max(T_ice(k)-T_pmp(k),0.0)

            ! Calculate basal mass balance as sum of all water produced in column,
            ! reset temperature to pmp  
            if (T_excess .gt. 0.0) then 
                melt_internal = melt_internal + T_excess * H_ice*(zeta_ac(k)-zeta_ac(k-1))*cp(k) / (L_ice * dt) 
                T_ice(k)      = T_pmp(k)
            end if 
            
        end do 

        ! Make sure base is below pmp too (mass/energy balance handled via bmb_grnd calculation externally)
        k = 1 
        if (T_ice(k) .gt. T_pmp(k)) T_ice(k) = T_pmp(k)


        ! Calculate basal mass balance (valid for grounded ice only)
        call calc_bmb_grounded(bmb_grnd,T_ice(1)-T_pmp(1),Q_ice_b,Q_b,Q_geo_now,f_grnd,rho_ice)
            
!         ! Get temperature gradient at ice base
!         if (H_ice .gt. 0.0_prec) then 
!             dz = H_ice * (zeta_aa(2)-zeta_aa(1))
!             dTdz_b = (T_ice(2) - T_ice(1)) / dz 
!         else 
!             dTdz_b = 0.0_prec 
!         end if 
        
!         ! Calculate basal mass balance (valid for grounded ice only)
!         call calc_bmb_grounded(bmb_grnd,T_ice(1)-T_pmp(1),dTdz_b,kt(1),rho_ice,Q_b,Q_geo_now,f_grnd)

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
        real(prec), intent(IN)    :: zeta_ac(:)    ! nz_ac == nz_aa-1 

        ! Local variables 
        integer :: k, nz_layers, nz_aa    

        nz_aa = size(zeta_aa)

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


