program test_icetemp 

    use ncio  
    
    use defs
    use thermodynamics 
    use icetemp 
    use iceage 

    implicit none 

    type icesheet_vectors
        real(prec), allocatable :: zeta(:)    ! [-] Sigma coordinates from 0:1 (either height or depth)
        real(prec), allocatable :: zeta_ac(:) ! nz-1 [-] sigma coordinates for internal ice layer edges
        real(prec), allocatable :: dzeta_a(:) ! nz [-] sigma coord helper for internal ice layer midpoints
        real(prec), allocatable :: dzeta_b(:) ! nz [-] sigma coord helper for internal ice layer midpoints
        real(prec), allocatable :: T_ice(:)   ! [K] Ice temperature 
        real(prec), allocatable :: T_pmp(:)   ! [K] Ice pressure melting point 
        real(prec), allocatable :: cp(:)      ! [] Ice heat capacity 
        real(prec), allocatable :: kt(:)      ! [] Ice conductivity  
        real(prec), allocatable :: uz(:)      ! [m a-1] Vertical velocity 
        real(prec), allocatable :: advecxy(:) ! [] Horizontal heat advection magnitude
        real(prec), allocatable :: Q_strn(:)  ! [K a-1] Strain heating 
        real(prec), allocatable :: t_dep(:)   ! [a] Deposition time 

    end type 

    type icesheet 
        real(prec) :: H_ice             ! [m] Ice thickness 
        real(prec) :: H_w               ! [m] Water present at the ice base 
        real(prec) :: T_srf             ! [K] Ice surface temperature 
        real(prec) :: smb               ! [m a**-1] Surface mass balance
        real(prec) :: bmb               ! [m a**-1] Basal mass balance
        real(prec) :: dTdz_b            ! [K m-1] Basal vertical temperature gradient 
        real(prec) :: Q_geo             ! [mW m-2] Geothermal heat flux 
        real(prec) :: Q_b               ! [] Basal heat production 
        logical    :: is_float          ! [-] Floating flag 
        character(len=56) :: age_method ! Method to use for age calculation 
        real(prec) :: age_impl_kappa    ! [m2 a-1] Artificial diffusion term for implicit age solving 
        type(icesheet_vectors) :: vec   ! For height coordinate systems with k=1 base and k=nz surface

    end type 

    ! Define different icesheet objects for use in proram
    type(icesheet) :: ice1
    type(icesheet) :: robin 
    type(icesheet) :: diff 
    
    ! Local variables
    real(prec)         :: t_start, t_end, dt, time  
    integer            :: n, ntot 
    character(len=512) :: file1D 
    real(prec)         :: dt_out 
    integer            :: nz
    logical            :: is_celcius 
    character(len=56)  :: age_method 
    real(prec)         :: age_impl_kappa
    real(prec)         :: t_base

    ! ===============================================================
    ! User options 

    t_start = 0.0       ! [yr]
    t_end   = 1000e3     ! [yr]
    dt      = 5.0       ! [yr]

    file1D  = "test.nc" 
    dt_out  = 10000.0      ! [yr] 

    nz      = 31           ! [--] Number of ice sheet points 

    is_celcius = .FALSE. 

    age_method     = "impl"     ! "expl" or "impl"
    age_impl_kappa = 1.5        ! [m2 a-1] Artificial diffusion for age tracing

    ! ===============================================================

    ! Initialize time and calculate number of time steps to iterate and 
    time = t_start 
    ntot = (t_end-t_start)/dt 
    
    ! Initialize icesheet object 
    call icesheet_allocate(ice1,nz=nz,zeta_scale="exp") 

    ! Prescribe initial eismint conditions for testing 
    call init_eismint_summit(ice1)
    !call init_k15_expa(ice1)

    ! Set age method and kappa 
    ice1%age_method     = age_method 
    ice1%age_impl_kappa = age_impl_kappa

    ! Calculate the robin solution for comparison 
    robin = ice1  
    robin%vec%T_ice = calc_temp_robin_column(robin%vec%zeta,robin%vec%T_pmp,robin%vec%kt,robin%vec%cp,rho_ice, &
                                       robin%H_ice,robin%T_srf,robin%smb,robin%Q_geo,robin%is_float)

    ! Write Robin solution 
    file1D = "robin.nc"
    call write_init(robin,filename=file1D,zeta=robin%vec%zeta,time_init=time)
    call write_step(robin,robin%vec,filename=file1D,time=time)

    ! Initialize output file for model and write intial conditions 
    file1D = "test.nc"
    call write_init(ice1,filename=file1D,zeta=ice1%vec%zeta,time_init=time)
    call write_step(ice1,ice1%vec,filename=file1D,time=time)

    ! Ensure zero basal water thickness to start 
    ice1%H_w = 0.0 

    ! Loop over time steps and perform thermodynamic calculations
    do n = 1, ntot 

        ! Get current time 
        time = t_start + n*dt 

        !if (time .ge. 100e3) ice1%T_srf = T0 - 5.0 
        !if (time .ge. 150e3) ice1%T_srf = T0 - 30.0 

        call calc_temp_column(ice1%vec%T_ice,ice1%bmb,ice1%dTdz_b,ice1%vec%T_pmp,ice1%vec%cp,ice1%vec%kt, &
                                ice1%vec%uz,ice1%vec%Q_strn,ice1%vec%advecxy,ice1%Q_b,ice1%Q_geo, &
                                ice1%T_srf,ice1%H_ice,ice1%H_w,ice1%is_float,ice1%vec%zeta, &
                                ice1%vec%zeta_ac,ice1%vec%dzeta_a,ice1%vec%dzeta_b,dt)

        ! Update basal water thickness 
        ice1%H_w = ice1%H_w - ice1%bmb*dt 

        ! Age calculations 
        t_base = ice1%vec%t_dep(1) 
        call calc_tracer_column(ice1%vec%t_dep,ice1%vec%uz,ice1%vec%advecxy*0.0,time,ice1%bmb, &
                                ice1%H_ice,ice1%vec%zeta,ice1%vec%zeta_ac,ice1%vec%dzeta_a,ice1%vec%dzeta_b, &
                                ice1%age_impl_kappa,dt)

        if (mod(time,dt_out)==0) then 
            call write_step(ice1,ice1%vec,filename=file1D,time=time,T_robin=robin%vec%T_ice)
        end if 

        if (mod(time,50.0)==0) then
            write(*,"(a,f14.4)") "time = ", time
        end if 

    end do

    write(*,*)
    write(*,*) "========================="
    write(*,*) 
    write(*,*) "Program finished."
    write(*,*)
    write(*,*) "========================="
    write(*,*)

contains 


    subroutine init_eismint_summit(ice)

        implicit none 

        type(icesheet), intent(INOUT) :: ice

        ! Local variables 
        integer :: k, nz, nz_ac   

        nz    = size(ice%vec%zeta)
        nz_ac = nz - 1 

        ! Assign point values
        ice%T_srf    = 239.0       ! [K] 
        ice%smb      = 0.1         ! [m/a]
        ice%bmb      = 0.0         ! [m/a]
        ice%Q_geo    = 42.0        ! [mW/m2]
        ice%H_ice    = 2997.0      ! [m] Summit thickness
        ice%H_w      = 0.0         ! [m] No basal water
        ice%Q_b      = 0.0         ! [] No basal frictional heating 
        ice%is_float = .FALSE.     ! Grounded point 

        ! EISMINT1
        ice%vec%cp      = 2009.0    ! [J kg-1 K-1]
        ice%vec%kt      = 6.67e7    ! [J a-1 m-1 K-1]
        
        ice%vec%Q_strn  = 0.0       ! [J a-1 m-3] No internal strain heating 
        ice%vec%advecxy = 0.0       ! [] No horizontal advection 

        ! Calculate pressure melting point 
        ice%vec%T_pmp = calc_T_pmp(ice%H_ice,ice%vec%zeta,T0) 

        if (is_celcius) then 
            ice%T_srf    = ice%T_srf    - T0
            ice%vec%T_pmp = ice%vec%T_pmp - T0 
        end if 

        ! Define initial temperature profile, linear 
        ! from T_srf at the surface to 10deg below freezing point at the base

        ice%vec%T_ice(nz) = ice%T_srf 
        ice%vec%T_ice(1)  = ice%vec%T_pmp(1) - 10.0 

        ! Intermediate layers are linearly interpolated 
        do k = 2, nz-1 
            ice%vec%T_ice(k) = ice%vec%T_ice(1)+ice%vec%zeta(k)*(ice%vec%T_ice(nz)-ice%vec%T_ice(1))
        end do 

        ! Define linear vertical velocity profile
        ice%vec%uz(nz_ac) = -ice%smb 
        ice%vec%uz(1)  = 0.0 
        do k = 2, nz_ac 
            ice%vec%uz(k) = ice%vec%uz(1)+(ice%vec%zeta_ac(k))*(ice%vec%uz(nz_ac)-ice%vec%uz(1))
        end do 
        
        return 

    end subroutine init_eismint_summit 

    subroutine init_k15_expa(ice)

        implicit none 

        type(icesheet), intent(INOUT) :: ice

        ! Local variables 
        integer :: k, nz 

        nz    = size(ice%vec%zeta)

        ! Assign point values
        ice%T_srf    = T0 - 30.0   ! [K]
        ice%smb      = 0.0         ! [m/a]
        ice%bmb      = 0.0         ! [m/a]
        ice%Q_geo    = 42.0        ! [mW/m2]
        ice%H_ice    = 1000.0      ! [m] Summit thickness
        ice%H_w      = 0.0         ! [m] No basal water
        ice%Q_b      = 0.0         ! [] No basal frictional heating 
        ice%is_float = .FALSE.     ! Grounded point 

        ! EISMINT1
        ice%vec%cp      = 2009.0    ! [J kg-1 K-1]
        ice%vec%kt      = 6.67e7    ! [J a-1 m-1 K-1]
        
        ice%vec%Q_strn  = 0.0       ! [] No internal strain heating 
        ice%vec%advecxy = 0.0       ! [] No horizontal advection 

        ! Calculate pressure melting point 
        ice%vec%T_pmp = calc_T_pmp(ice%H_ice,ice%vec%zeta,T0) 

        if (is_celcius) then 
            ice%T_srf    = ice%T_srf      - T0
            ice%vec%T_pmp = ice%vec%T_pmp - T0 
        end if 

        ! Define initial temperature profile
        ! (constant equal to surface temp)
        ice%vec%T_ice = ice%T_srf 
        where(ice%vec%T_ice .gt. T0) ice%vec%T_ice = T0 

        ! Intermediate layers are linearly interpolated 
        do k = 2, nz-1 
            ice%vec%T_ice(k) = ice%vec%T_ice(1)+ice%vec%zeta(k)*(ice%vec%T_ice(nz)-ice%vec%T_ice(1))
        end do 

        ! Define linear vertical velocity profile
        ice%vec%uz(nz) = -ice%smb 
        ice%vec%uz(1)  = 0.0 
        do k = 2, nz-1 
            ice%vec%uz(k) = ice%vec%uz(1)+(ice%vec%zeta(k))*(ice%vec%uz(nz)-ice%vec%uz(1))
        end do 
        
        return 

    end subroutine init_k15_expa 
    
    subroutine icesheet_allocate(ice,nz,zeta_scale)
        ! Allocate the ice sheet object 

        implicit none 

        type(icesheet), intent(INOUT) :: ice 
        integer,        intent(IN)    :: nz     ! Number of ice points (aa-nodes)
        character(*),   intent(IN)    :: zeta_scale

        ! Local variables 
        integer :: k, nz_ac 

        nz_ac = nz-1 

        ! First allocate 'up' variables (with vertical coordinate as height)

        ! Make sure all vectors are deallocated
        if (allocated(ice%vec%zeta))     deallocate(ice%vec%zeta)
        if (allocated(ice%vec%zeta_ac))  deallocate(ice%vec%zeta_ac)
        if (allocated(ice%vec%dzeta_a))  deallocate(ice%vec%dzeta_a)
        if (allocated(ice%vec%dzeta_b))  deallocate(ice%vec%dzeta_b)
        
        if (allocated(ice%vec%T_ice))   deallocate(ice%vec%T_ice)
        if (allocated(ice%vec%T_pmp))   deallocate(ice%vec%T_pmp)
        if (allocated(ice%vec%cp))      deallocate(ice%vec%cp)
        if (allocated(ice%vec%kt))      deallocate(ice%vec%kt)
        if (allocated(ice%vec%uz))      deallocate(ice%vec%uz)
        if (allocated(ice%vec%advecxy)) deallocate(ice%vec%advecxy)
        if (allocated(ice%vec%Q_strn))  deallocate(ice%vec%Q_strn)

        if (allocated(ice%vec%t_dep))   deallocate(ice%vec%t_dep)

        ! Allocate vectors with desired lengths
        allocate(ice%vec%zeta(nz))
        allocate(ice%vec%zeta_ac(nz_ac))
        allocate(ice%vec%dzeta_a(nz))
        allocate(ice%vec%dzeta_b(nz))
        
        allocate(ice%vec%T_ice(nz))
        allocate(ice%vec%T_pmp(nz))
        allocate(ice%vec%cp(nz))
        allocate(ice%vec%kt(nz))
        allocate(ice%vec%advecxy(nz))
        allocate(ice%vec%Q_strn(nz))

        allocate(ice%vec%uz(nz_ac))
        
        allocate(ice%vec%t_dep(nz))

        ! Initialize zeta 
!         ice%vec%zeta = 0.0  
!         do k = 1, nz 
!             ice%vec%zeta(k) = real(k-1,prec) / real(nz-1,prec)
!             !write(*,*) ice%vec%zeta(k)
!         end do  

        ice%vec%zeta_ac = 0.0  
        do k = 1, nz_ac 
            ice%vec%zeta_ac(k) = real(k-1,prec) / real(nz_ac-1,prec)
            !write(*,*) ice%vec%zeta(k)
        end do 

        ! Scale zeta to produce different resolution through column if desired
        ! zeta_scale = ["linear","exp","wave"]
        select case(trim(zeta_scale))
            
            case("exp")
                ! Increase resolution at the base 
                ice%vec%zeta_ac = ice%vec%zeta_ac**2.0

            case("tanh")
                ! Increase resolution at base and surface 

                ice%vec%zeta_ac = tanh(1.5*pi*(ice%vec%zeta_ac-0.65))
                ice%vec%zeta_ac = ice%vec%zeta_ac - minval(ice%vec%zeta_ac)
                ice%vec%zeta_ac = ice%vec%zeta_ac / maxval(ice%vec%zeta_ac)

            case DEFAULT
            ! Do nothing, scale should be linear as defined above
        
        end select  

        ! Calculate zeta_ac (zeta on ac-nodes)
        !ice%vec%zeta_ac = calc_zeta_ac(ice%vec%zeta)

        ice%vec%zeta(1) = 0.0 
        do k = 2, nz-1
            ice%vec%zeta(k) = 0.5 * (ice%vec%zeta_ac(k-1) + ice%vec%zeta_ac(k))
        end do 
        ice%vec%zeta(nz) = 1.0

!         ! =======================================================
!         ice%vec%zeta_ac = 0.0  
!         do k = 1, nz_ac 
!             ice%vec%zeta_ac(k) = real(k-1,prec) / real(nz_ac-1,prec)
!             !write(*,*) ice%vec%zeta(k)
!         end do  

!         ice%vec%zeta(1) = 0.0 
!         do k = 2, nz-1
!             ice%vec%zeta(k) = 0.5 * (ice%vec%zeta_ac(k-1) + ice%vec%zeta_ac(k))
!         end do 
!         ice%vec%zeta(nz) = 1.0 
!         ! =======================================================

        ! Define thermodynamic zeta helper derivative variables dzeta_a/dzeta_b
        call calc_dzeta_terms(ice%vec%dzeta_a,ice%vec%dzeta_b,ice%vec%zeta,ice%vec%zeta_ac)

!         do k = 1, nz 
!             write(*,*) ice%vec%zeta(k), ice%vec%dzeta_a(k), ice%vec%dzeta_b(k) 
!         end do 

!         write(*,*) "-----"

!         do k = 1, nz-1 
!             write(*,*) ice%vec%zeta_ac(k)
!         end do 

!         stop 
        
        ! Initialize remaining vectors to zero 
        ice%vec%T_ice   = 0.0 
        ice%vec%T_pmp   = 0.0 
        ice%vec%cp      = 0.0
        ice%vec%kt      = 0.0
        ice%vec%uz      = 0.0
        ice%vec%advecxy = 0.0  
        ice%vec%Q_strn  = 0.0

        ice%vec%t_dep   = 0.0 

        write(*,*) "Allocated icesheet variables."

        return 

    end subroutine icesheet_allocate 

    subroutine write_init(ice,filename,zeta,time_init)

        implicit none 

        type(icesheet),   intent(IN) :: ice 
        character(len=*), intent(IN) :: filename 
        real(prec),       intent(IN) :: zeta(:)  
        real(prec),       intent(IN) :: time_init

        ! Initialize netcdf file and dimensions
        call nc_create(filename)
        call nc_write_dim(filename,"zeta",  x=zeta,    units="1")
        call nc_write_dim(filename,"time",  x=time_init,dx=1.0_prec,nx=1,units="years",unlimited=.TRUE.)

        return

    end subroutine write_init 
    
    subroutine write_step(ice,vecs,filename,time,T_robin)

        implicit none 
        
        type(icesheet),         intent(IN) :: ice
        type(icesheet_vectors), intent(IN) :: vecs
        character(len=*),       intent(IN) :: filename
        real(prec),             intent(IN) :: time
        real(prec), optional,   intent(IN) :: T_robin(:) 

        ! Local variables
        integer    :: ncid, n
        real(prec) :: time_prev 
        character(len=12), parameter :: vert_dim = "zeta"

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
        if (abs(time-time_prev).gt.1e-5) n = n+1 

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! Update variables (vectors) 
        call nc_write(filename,"T_ice",  vecs%T_ice,  units="K",   long_name="Ice temperature",           dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"T_pmp",  vecs%T_pmp,  units="",    long_name="Ice pressure melting point",dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"T_prime",vecs%T_ice-vecs%T_pmp,units="K",long_name="Ice temperature",     dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        
        call nc_write(filename,"cp",     vecs%cp,     units="J kg-1 K-1",   long_name="Ice heat capacity",       dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"kt",     vecs%kt,     units="J a-1 m-1 K-1",long_name="Ice thermal conductivity",dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"uz",     vecs%uz,     units="m a**-1",  long_name="Ice vertical velocity",   dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"advecxy",vecs%advecxy,units="",         long_name="Ice horizontal advection",dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"Q_strn", vecs%Q_strn, units="J a-1 m-3",long_name="Ice strain heating",      dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        
        call nc_write(filename,"t_dep", vecs%t_dep, units="a",long_name="Deposition time",      dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        
        ! Update variables (points) 
        call nc_write(filename,"smb",     ice%smb,units="m a**-1",long_name="Surface mass balance",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"bmb",     ice%bmb,units="m a**-1",long_name="Basal mass balance",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"Q_b",     ice%Q_b,units="",long_name="Basal heating",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"Q_geo",   ice%Q_geo,units="mW m**-2",long_name="Geothermal heat flux",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"T_srf",   ice%T_srf,units="K",long_name="Surface temperature",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"H_ice",   ice%H_ice,units="m",long_name="Ice thickness",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"H_w",     ice%H_w,units="m",long_name="Basal water thickness",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"is_float",ice%is_float,units="",long_name="Floating flag",dim1="time",start=[n],ncid=ncid)
        
        ! If available, compare with Robin analytical solution 
        if (present(T_robin)) then 
            call nc_write(filename,"T_diff",  vecs%T_ice-T_robin,  units="K", &
                            long_name="Ice temperature difference with Robin solution", &
                            dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        end if 

        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step

    function calc_zeta_ac(zeta_aa) result(zeta_ac)
        ! Calculate the vertical layer-edge axis (vertical ac-nodes)
        ! given the vertical layer-boundary axis as input 

        implicit none 

        real(prec), intent(IN)  :: zeta_aa(:) 
        real(prec) :: zeta_ac(size(zeta_aa))

        ! Local variables
        integer :: k, nz_aa, nz_ac 

        nz_aa  = size(zeta_aa)
        nz_ac  = nz_aa-1 

        ! Get zeta_ac (edges of zeta layers - between zeta_aa values)
        zeta_ac(1) = 0.0 
        do k = 2, nz_ac-1
            zeta_ac(k) = 0.5 * (zeta_aa(k)+zeta_aa(k+1))
        end do 
        zeta_ac(nz_ac) = 1.0 
        
        if (.FALSE.) then 
            ! zeta_ac(nz_ac==nz_aa) solution for alternate bottom boundary condition

            ! Get zeta_ac (edges of zeta layers - between zeta_aa values) 
            do k = 1, nz_ac-1
                zeta_ac(k) = 0.5 * (zeta_aa(k)+zeta_aa(k+1))
            end do 
            zeta_ac(nz_ac) = 1.0 

        end if 

        return 

    end function calc_zeta_ac

end program test_icetemp 
