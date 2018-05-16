module icetemp_grisli
    ! Wrapping the original grisli icetemp solution for yelmo 

    use defs, only : prec, pi, g, sec_year, T0, rho_ice, rho_sw, rho_w
    use solver_tridiagonal, only : tridiag 
    use thermodynamics 

    implicit none
    
    ! Impose parameter choice here for testing conductive bedrock or not
    logical,    parameter :: conductive_bedrock = .FALSE. 
    
    private
    public :: calc_icetemp_grisli_column

contains 
    
    subroutine calc_icetemp_grisli_column(ibase,T_ice,T_rock,T_pmp,cp,ct,uz,Q_strn,advecxy,Q_b, &
                                            Q_geo,T_srf,H_ice,H_w,is_float,dt)
        ! GRISLI solver for thermodynamics for a given column of ice 

        ! Note T_ice, T_pmp in [degC], not [K]

        implicit none 

        integer,    intent(INOUT) :: ibase
        real(prec), intent(OUT) :: T_ice(:)
        real(prec), intent(OUT) :: T_rock(:)
        real(prec), intent(OUT) :: T_pmp(:)
        real(prec), intent(OUT) :: cp(:)
        real(prec), intent(OUT) :: ct(:)
        real(prec), intent(IN)  :: uz(:)
        real(prec), intent(IN)  :: Q_strn(:)
        real(prec), intent(IN)  :: advecxy(:)
        real(prec), intent(IN)  :: Q_b 
        real(prec), intent(IN)  :: Q_geo 
        real(prec), intent(IN)  :: T_srf
        real(prec), intent(IN)  :: H_ice 
        real(prec), intent(IN)  :: H_w
        logical,    intent(IN)  :: is_float
        real(prec), intent(IN)  :: dt 

        ! Local variables 
        integer    :: k, nz   
        integer    :: nzm, nzz, nfracq  
        real(prec) :: dzm, de, da, ro, rom, romg, row, cl, cm, cpm
        real(prec) :: dzz, dah, dou, dzi  
        real(prec) :: ct_bas, ct_haut 
        real(prec) :: ctm, fracq 
        real(prec) :: acof1, bcof1, ccof1, s0mer, tbmer 
        real(prec) :: tbdot, tdot, tss
        real(prec) :: bmelt
        real(prec) :: Q_geo_now   
        integer    :: ifail 

        real(prec), allocatable :: aa(:), bb(:), cc(:), rr(:), hh(:) 
        real(prec), allocatable :: abis(:), bbis(:), cbis(:), rbis(:), hbis(:)  
        real(prec), allocatable :: T(:), T_new(:) 

        nz = size(T_ice,1) 

        nzm     = size(T_rock,1)   ! Number of bedrock points 
        nzz     = nz+nzm           ! Total grid points (ice plus bedrock)
        dzm     = 600.0            ! [m] Bedrock step height 
        de      = 1.0/(nz-1)       ! vertical step in ice and mantle
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
        bcof1   = 0.0901
        ccof1   = 7.61e-4
        s0mer   = 34.75

        ! Convert units of Geothermal heat flux
        Q_geo_now = -Q_geo *sec_year*1e-3   ! [mW/m2] => [J a-1 m-2]

        ! Derived constants
        ctm     = dt*cm/rom/cpm/dzm/dzm  

        nfracq = 1
        fracq  = (1.0-(1.-de/2.0)**nfracq)/nfracq

        ! Allocate local variables 
        allocate(aa(nzz),bb(nzz),cc(nzz),rr(nzz),hh(nzz))
        allocate(abis(nzz),bbis(nzz),cbis(nzz),rbis(nzz),hbis(nzz))

        allocate(T(nzz),T_new(nzz))

        ! Populate local temperature column (ice+rock)
        T(1:nz)        = T_ice 
        T(nz+1:nz+nzm) = T_rock 

        ! Populate new solution to be safe 
        T_new = T 

        ! GRISLI THERMAL PROPERTIES ===========================================
!         ! Calculate thermal properties
!         do k = 1, nz  
!             T_pmp(k) = -0.00087*(k-1)*de*H_ice               ! [degC]
!             cp(k)    = (2115.3+7.79293*T_ice(k))*ro          ! [J m-3 K-1]
!             ct(k)    = 3.1014e8*exp(-0.0057*(T_ice(k)+T0))   ! [J m-1 K-1 a-1] 
!         end do 
        ! =====================================================================

        ! EISMINT THERMAL PROPERTIES ==========================================
        cp = 2009.0 * ro   ! [J m-3 K-1]
        ct = 6.6e7         ! [J m-1 K-1 a-1] 
        ! =====================================================================
        
        ! Calculate the sea temperature below ice shelf if point is floating 
        if (is_float) then
            tbmer = acof1*s0mer + bcof1 + ccof1*H_ice*ro/row
        else 
            tbmer = 0.0 
        end if 

        ! Diagnose basal melt rate and basal state (temperate, frozen, etc.)
        call bmelt_grounded_column(bmelt,ibase,T,ct,Q_b,H_ice,H_w,Q_geo_now,is_float, &
                                    de,dzm,cm,ro,cl,dt)


        ! In the bedrock, the elements of the tridiagonal matrix are invariant, 
        ! so populate them now 
        do k=nz+1,nz+nzm-1
            aa(k) = -ctm
            bb(k) = 1.0+2.0*ctm
            cc(k) = -ctm
        end do

        ! Boundary conditions at the base of the bedrock 
        aa(nz+nzm) = -1.0
        bb(nz+nzm) =  1.0
        cc(nz+nzm) =  0.0
        rr(nz+nzm) = -dzm*Q_geo_now/cm


        if (H_ice .gt. 10.0) then
            ! General case (H_ice>10m)

            ! =============================================================
            ! Ice sheet conditions

            ! Prescribe surface temperature as boundary condition
            T(1) = min(0.0,T_srf)

            aa(1) = 0.0
            bb(1) = 1.0
            cc(1) = 0.0
            rr(1) = T(1)

            ! Internal ice sheet points
            dou    = dt/de/de/H_ice/H_ice
            dah    = dt/H_ice/de
            ct_bas = 2.0*(Ct(1)*Ct(2))/(Ct(1)+Ct(2))

            do k = 2, nz-1

                dzz = dou/cp(k)

                ! Conductivity as the Harmonic average of the grid points
                Ct_haut = ct_bas
                Ct_bas  = 2.0*(Ct(k)*Ct(k+1))/(Ct(k)+Ct(k+1))

                ! Vertical advection (centered)
                aa(k) = -dzz*ct_haut-uz(k)*dah/2.0
                bb(k) = 1.0+dzz*(ct_haut+ct_bas)
                cc(k) = -dzz*ct_bas+uz(k)*dah/2.0
                rr(k) = T(k)+dt*Q_strn(k)-dt*advecxy(k)

            end do
            

            ! Limits at the base of the ice sheet 

            ! Frozen base 
            if (.not. is_float .and. ((ibase.eq.1).or.(ibase.eq.4) &
                .or. ((ibase.eq.5).and.(T(nz).lt.T_pmp(nz)))) ) then

                 ibase  = 1
                 bb(nz) = 1.0

                if (conductive_bedrock) then
                    ! With conductive bedrock 
                    dzi    = H_ice*de*cm/ct(nz)
                    aa(nz) = -dzm/(dzm+dzi)
                    cc(nz) = -dzi/(dzm+dzi)
                    rr(nz) = H_ice*de*Q_b*dzm/ct(nz)/(dzm+dzi)

                else
                    ! Without conductive bedrock 
                    aa(nz) = -1.0
                    cc(nz) =  0.0
                    rr(nz) = -(Q_geo_now-Q_b)/ct(nz)*H_ice*de

               endif

            else
                ! Temperate ice or shelf ice 

                 if (ibase.eq.5) ibase = 2
                 ibase = max(ibase,2)
                 aa(nz) = 0.0
                 bb(nz) = 1.0
                 cc(nz) = 0.0

                if (.not. is_float) then
                    rr(nz) = T_pmp(nz)
                else
                    rr(nz) = Tbmer
                end if

            endif

            ! =============================================================
            ! Bedrock conditions and solver

            if (conductive_bedrock) then
                ! With conductive bedrock 

                do k=nz+1,nz+nzm-1
                    rr(k)=T(k)
                end do
                
                call tridiag(aa,bb,cc,rr,hh,nz+nzm,ifail)

            else
                ! Without conductive bedrock 

                call tridiag(aa,bb,cc,rr,hh,nz,ifail)

                ! Prescribe solution for bedrock points (linear with ghf)
                do k=nz+1,nz+nzm
                    hh(k) = hh(nz)-dzm*(k-nz)*Q_geo_now/cm
                end do

            endif

            ! Check tridiag results
            if (ifail.eq.1) then
                write(*,*) 'Error in tridiag solver.'
                write(*,*) 'a, b, c, r, u'
                do k = 1, nz+nzm
                    write(*,*)'k = ', k
                    write(*,*) aa(k),bb(k),cc(k),rr(k),hh(k)
                end do
                stop
            end if

            ! Fill in new temperature solution 
            do k=1,nz+nzm
              T_new(k) = hh(k)
            end do
            
            ! Diagnose rate of temperature change at ice base 
            tbdot = (T_new(nz)-T(nz))/dt


        else if (H_ice .le. 10.0) then
            ! For very thin ice or no ice:
            ! To avoid problems with very thin ice layers, prescribe
            ! a linear profile according to the geothermal heat flux
            ! Points with no ice are treated the same way 

            if (conductive_bedrock .and. .not. is_float) then
                ! With conductive bedrock 

                if (H_ice .gt. 0.0) then
                    ! Ice exists, impose the gradient from the bedrock 

                    dou = (T(nz+1)-T(nz))/dzm*cm
                    dou = dou/ct(nz)*de*H_ice

                    tss = min(0.,T_srf)
                    do k = 1, nz
                        T_new(k) = tss+dou*(k-1.0)
                    end do

                else
                    ! Ice-free point 

                    tss = T_srf
                    do k = 1, nz 
                        T_new(k) = tss
                    end do
                end if

                ! Perform calculations in the bedrock even if no ice exists

                ! Ice base (bedrock surface?)
                aa(nz) = 0.0 
                bb(nz) = 1.0 
                cc(nz) = 0.0 
                rr(nz) = T_new(nz)

                ! Internal bedrock layers 
                do k = nz+1, nz+nzm-1
                    aa(k) = -ctm
                    bb(k) = 1.0+2.0*ctm
                    cc(k) = -ctm
                    rr(k) = T(k)
                end do

                ! Populate table just for solving bedrock     
                do k = 1, nzm+1
                    abis(k)=aa(nz-1+k)
                    bbis(k)=bb(nz-1+k)
                    cbis(k)=cc(nz-1+k)
                    rbis(k)=rr(nz-1+k)
                end do

                ! Solve bedrock 
                call tridiag(abis,bbis,cbis,rbis,hbis,nzm+1,ifail)

                ! Check tridiag solver results
                IF (ifail .eq. 1) then 
                    write(*,*)'Error in tridiag solver: bedrock only.'
                    write(*,*) 'a, b, c, r, u'
                    do k=1,nz+nzm 
                        write(*,*)'k = ',k
                        write(*,*) abis(k),bbis(k),cbis(k),rbis(k),hbis(k)
                    end do
                    stop 

                end if 

                ! Repopulate solution for bedrock points
                do k=1, nzm+1
                    hh(nz-1+k) = hbis(k)
                end do

            else
                ! Without conductive bedrock 

                                         
                if (H_ice .gt. 0.0 .and. .not. is_float) then
                    ! Grounded ice sheet point

                    dou = -Q_geo_now/ct(nz)*de*H_ice
                    tss = min(0.0,T_srf)
                    do k = 1, nz 
                        T_new(k) = tss+dou*(k-1.0)
                    end do

                else if (H_ice .gt. 0.0 .and. is_float) then
                    ! Shelf ice sheet point 

                    tss = min(0.0,T_srf)
                    dou = (tbmer-tss)*de
                    do k = 1, nz 
                        T_new(k) = tss+dou*(k-1.0)
                    end do

                else
                    ! Ice-free point 

                    tss = T_srf 
                    do k = 1, nz 
                        T_new(k) = tss
                    end do

                end if

                ! Calculate temperature in the bedrock as a linear gradient of ghf
                if (.not. is_float) then 
                    ! Grounded point 

                    do k = nz+1, nz+nzm 
                        hh(k) = T_new(nz)-dzm*(k-nz)*Q_geo_now/cm
                    end do

                else
                    ! Floating point 

                    do k = nz+1, nz+nzm 
                        hh(K) = Tbmer-dzm*(k-nz)*Q_geo_now/cm
                    end do
                
                end if 


            end if
            

            ! Populate new solution in bedrock 
            do k = nz+1, nz+nzm 
                T_new(k) = hh(k)
            end do

            ! Diagnose rate of change of temperature at ice base
            tbdot = (T_new(nz) - T(nz)) / dt 

            bmelt = 0.0
            ibase = 5
            !Q_b   = 0.0     ! ajr: now calculated externally

        end if

        ! Ensure T_pmp at the ice surface is equal to zero [degC]
        T_pmp(1) = 0.0


        ! Apply limits: less than 3degC / 5 yrs variation and no colder than -70 degC
        do k=1,nz+nzm

            if (H_ice .gt. 10.0) then

                tdot = (T_new(k)-T(k))/dt

                if (k .gt. 1 .and. tdot .lt. -3.0) tdot = -3.0
                if (k .gt. 1 .and. tdot .gt.  3.0) tdot =  3.0

                T(k) = T(k)+tdot

                if (T(k) .lt. -70.0) T(k) = -70.0

            end if

        end do

        ! Limit internal ice temperatures to the pressure melting point 
        do k = 2, nz

            if (T(k) .gt. T_pmp(k)) then
                T(k)  = T_pmp(k)
                ibase = 2
            end if

        end do

        ! Ensure ice temperature for very thin/ no ice points is equal
        ! to surface temperature 
        if (H_ice .le. 1.0) T(1:nz) = T_srf 

        ! Store solution back in output variables 
        T_ice  = T(1:nz)
        T_rock = T(nz+1:nz+nzm)

        return 

    end subroutine calc_icetemp_grisli_column 


    subroutine bmelt_grounded_column(bmelt,ibase,T,ct,Q_b,H_ice,H_w,Q_geo_now,is_float,de,dzm,cm,ro,cl,dt)
        ! Diagnose the basal melting rate and the state of the basal 
        ! ice (temperate, frozen, etc), for internal use in calc_icetemp_grisli_column
        
        implicit none 

        real(prec), intent(OUT)   :: bmelt 
        integer,    intent(INOUT) :: ibase 
        real(prec), intent(IN)    :: T(:) 
        real(prec), intent(IN)    :: ct(:) 
        real(prec), intent(IN)    :: Q_b
        real(prec), intent(IN)    :: H_ice 
        real(prec), intent(IN)    :: H_w 
        real(prec), intent(IN)    :: Q_geo_now    ! [J a-1 m-2]
        logical,    intent(IN)    :: is_float 
        real(prec), intent(IN)    :: de
        real(prec), intent(IN)    :: dzm 
        real(prec), intent(IN)    :: cm 
        real(prec), intent(IN)    :: ro 
        real(prec), intent(IN)    :: cl 
        real(prec), intent(IN)    :: dt 

        ! Local variables 
        integer :: nz 

        nz = size(ct) 

        if (.not. is_float .and. H_ice .gt. 10.0 .and. ibase .ne. 1) then 

            if (conductive_bedrock) then
                ! With conductive bedrock

                bmelt = (ct(nz)*(T(nz-1)-T(nz))/de/H_ice &
                        - cm*(T(nz)-T(nz+1))/dzm+Q_b)/ro/cl

            else
                ! Without conductive bedrock 

                bmelt = (ct(nz)*(T(nz-1)-T(nz))/de/H_ice &
                        + (Q_b-Q_geo_now))/ro/cl
            
            end if


            if (bmelt .ge. 0.0) then
                ! Basal melt present, temperate base
                
                ibase = 2
                                   
            else
                ! Refreezing 

                if (H_w .gt. -bmelt*dt) then
                    ! If basal water will remain after refreezing, still temperate 
                    
                    ibase=2

                else if (ibase .eq. 3) then
                    ! Point became frozen now  
                    ibase = 4
                    bmelt = 0.0
                    !H_w   = 0.0    ! ajr: now handled externally 

                else if (H_w .le. -bmelt*dt) then
                    ! Refreezing, but there is not enough water

                    ibase = 3
                
                end if

            end if

        else
            ! Thin or no ice 

            bmelt = 0.0

        end if 

        return

    end subroutine bmelt_grounded_column 

end module icetemp_grisli 

