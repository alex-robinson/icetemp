program test_icetemp 

    use ncio  
    
    use defs
    use thermodynamics 
    use icetemp 

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
        real(prec), allocatable :: uz(:)      ! [] Vertical velocity 
        real(prec), allocatable :: advecxy(:) ! [] Horizontal heat advection magnitude
        real(prec), allocatable :: Q_strn(:)  ! [] Strain heating 
        
    end type 

    type icesheet 
        real(prec) :: H_ice     ! [m] Ice thickness 
        real(prec) :: H_w       ! [m] Water present at the ice base 
        real(prec) :: T_srf     ! [K] Ice surface temperature 
        real(prec) :: smb       ! [m a**-1] Surface mass balance
        real(prec) :: bmb       ! [m a**-1] Basal mass balance
        real(prec) :: dTdz_b    ! [K m-1] Basal vertical temperature gradient 
        real(prec) :: Q_geo     ! [mW m-2] Geothermal heat flux 
        real(prec) :: Q_b       ! [] Basal heat production 
        logical    :: is_float  ! [-] Floating flag 
        
        type(icesheet_vectors) :: up     ! For height coordinate systems with k=1 base and k=nz surface

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
    integer            :: nz, nzt
    logical            :: is_celcius 

    ! ===============================================================
    ! User options 

    t_start = 0.0       ! [yr]
    t_end   = 200e3     ! [yr]
    dt      = 5.0       ! [yr]

    file1D  = "test.nc" 
    dt_out  = 10000.0      ! [yr] 

    nz      = 21           ! [--] Number of ice sheet points 

    is_celcius = .FALSE. 
    ! ===============================================================

    ! Initialize time and calculate number of time steps to iterate and 
    time = t_start 
    ntot = (t_end-t_start)/dt 
    
    ! Initialize icesheet object 
    call icesheet_allocate(ice1,nz=nz)
    nzt = nz - 1 

    ! Prescribe initial eismint conditions for testing 
    call init_eismint_summit(ice1)

    ! Initialize output file and write intial conditions 
    call write_init(ice1,filename=file1D,zeta=ice1%up%zeta,time_init=time)
    call write_step(ice1,ice1%up,filename=file1D,time=time)
    
    ! Loop over time steps and perform thermodynamic calculations
    do n = 1, ntot 

        ! Get current time 
        time = t_start + n*dt 

        call calc_temp_column(ice1%up%T_ice,ice1%bmb,ice1%dTdz_b,ice1%up%T_pmp,ice1%up%cp,ice1%up%kt, &
                                ice1%up%uz,ice1%up%Q_strn,ice1%up%advecxy,ice1%Q_b,ice1%Q_geo, &
                                ice1%T_srf,ice1%H_ice,ice1%H_w,ice1%is_float,ice1%up%zeta, &
                                ice1%up%zeta_ac,ice1%up%dzeta_a,ice1%up%dzeta_b,dt)

        if (mod(time,dt_out)==0) then 
            call write_step(ice1,ice1%up,filename=file1D,time=time)
        end if 

        if (mod(time,50.0)==0) then
            write(*,"(a,f14.4)") "time = ", time
        end if 

    end do 

    ! Also calculate the robin solution for comparison 
    robin = ice1  
    robin%up%T_ice = calc_temp_robin_column(robin%up%zeta,robin%up%T_pmp,robin%up%kt,robin%up%cp,rho_ice, &
                                       robin%H_ice,robin%T_srf,robin%smb,robin%Q_geo,robin%is_float)

    ! Write Robin solution for comparison 
    file1D = "robin.nc"
    call write_init(robin,filename=file1D,zeta=robin%up%zeta,time_init=time)
    call write_step(robin,robin%up,filename=file1D,time=time)

    ! Compare our solution with robin and write comparison results 
    diff = robin 
    diff%up%T_ice    = ice1%up%T_ice - robin%up%T_ice  

    file1D = "diff.nc"
    call write_init(diff,filename=file1D,zeta=diff%up%zeta,time_init=time)
    call write_step(diff,diff%up,filename=file1D,time=time)


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

        nz    = size(ice%up%zeta)
        nz_ac = size(ice%up%zeta_ac) 

        ! Assign point values
        ice%T_srf    = 239.0       ! [K] 
        ice%smb      = 0.5         ! [m/a]
        ice%bmb      = 0.0         ! [m/a]
        ice%Q_geo    = 42.0        ! [mW/m2]
        ice%H_ice    = 2997.0      ! [m] Summit thickness
        ice%H_w      = 0.0         ! [m] No basal water
        ice%Q_b      = 0.0         ! [] No basal frictional heating 
        ice%is_float = .FALSE.     ! Grounded point 

        ! EISMINT1
        ice%up%cp      = 2009.0    ! [J kg-1 K-1]
        ice%up%kt      = 6.67e7    ! [J a-1 m-1 K-1]
        
        ice%up%Q_strn  = 0.0       ! [] No internal strain heating 
        ice%up%advecxy = 0.0       ! [] No horizontal advection 

        ! Calculate pressure melting point 
        ice%up%T_pmp = calc_T_pmp(ice%H_ice,ice%up%zeta,T0) 

        if (is_celcius) then 
            ice%T_srf    = ice%T_srf    - T0
            ice%up%T_pmp = ice%up%T_pmp - T0 
        end if 

        ! Define initial temperature profile, linear 
        ! from T_srf at the surface to 10deg below freezing point at the base

        ice%up%T_ice(nz) = ice%T_srf 
        ice%up%T_ice(1)  = ice%up%T_pmp(1) - 10.0 

        ! Intermediate layers are linearly interpolated 
        do k = 2, nz-1 
            ice%up%T_ice(k) = ice%up%T_ice(1)+ice%up%zeta(k)*(ice%up%T_ice(nz)-ice%up%T_ice(1))
        end do 

        ! Define linear vertical velocity profile
        ice%up%uz(nz) = -ice%smb 
        ice%up%uz(1)  = 0.0 
        do k = 2, nz-1 
            ice%up%uz(k) = ice%up%uz(1)+(ice%up%zeta(k))*(ice%up%uz(nz)-ice%up%uz(1))
        end do 
        
        return 

    end subroutine init_eismint_summit 

    subroutine icesheet_allocate(ice,nz)
        ! Allocate the ice sheet object 

        implicit none 

        type(icesheet), intent(INOUT) :: ice 
        integer,        intent(IN)    :: nz     ! Number of ice points (aa-nodes)

        ! Local variables 
        integer :: k, nz_ac 

        nz_ac = nz+1 

        ! First allocate 'up' variables (with vertical coordinate as height)

        ! Make sure all vectors are deallocated
        if (allocated(ice%up%zeta))     deallocate(ice%up%zeta)
        if (allocated(ice%up%zeta_ac))  deallocate(ice%up%zeta_ac)
        if (allocated(ice%up%dzeta_a))  deallocate(ice%up%dzeta_a)
        if (allocated(ice%up%dzeta_b))  deallocate(ice%up%dzeta_b)
        
        if (allocated(ice%up%T_ice))   deallocate(ice%up%T_ice)
        if (allocated(ice%up%T_pmp))   deallocate(ice%up%T_pmp)
        if (allocated(ice%up%cp))      deallocate(ice%up%cp)
        if (allocated(ice%up%kt))      deallocate(ice%up%kt)
        if (allocated(ice%up%uz))      deallocate(ice%up%uz)
        if (allocated(ice%up%advecxy)) deallocate(ice%up%advecxy)
        if (allocated(ice%up%Q_strn))  deallocate(ice%up%Q_strn)

        ! Allocate vectors with desired lengths
        allocate(ice%up%zeta(nz))
        allocate(ice%up%zeta_ac(nz_ac))
        allocate(ice%up%dzeta_a(nz))
        allocate(ice%up%dzeta_b(nz))
        
        allocate(ice%up%T_ice(nz))
        allocate(ice%up%T_pmp(nz))
        allocate(ice%up%cp(nz))
        allocate(ice%up%kt(nz))
        allocate(ice%up%uz(nz))
        allocate(ice%up%advecxy(nz))
        allocate(ice%up%Q_strn(nz))

        ! Initialize zeta 
        ice%up%zeta = 0.0  
        do k = 1, nz 
            ice%up%zeta(k) = real(k-1,prec) / real(nz-1,prec)
            !write(*,*) ice%up%zeta(k)
        end do 

        ! Nonlinear zeta (smaller dz increments at the base):
        ice%up%zeta = ice%up%zeta**2.0
        
        ! Calculate zeta_ac (zeta on ac-nodes)
        ice%up%zeta_ac = calc_zeta_ac(ice%up%zeta)

        ! Define thermodynamic zeta helper derivative variables dzeta_a/dzeta_b
        call calc_dzeta_terms(ice%up%dzeta_a,ice%up%dzeta_b,ice%up%zeta,ice%up%zeta_ac)

        ! Initialize remaining vectors to zero 
        ice%up%T_ice   = 0.0 
        ice%up%T_pmp   = 0.0 
        ice%up%cp      = 0.0
        ice%up%kt      = 0.0
        ice%up%uz      = 0.0
        ice%up%advecxy = 0.0  
        ice%up%Q_strn  = 0.0

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
        call nc_write_dim(filename,"zeta",    x=zeta,    units="1")
        call nc_write_dim(filename,"time",  x=time_init,dx=1.0_prec,nx=1,units="years",unlimited=.TRUE.)

        return

    end subroutine write_init 
    
    subroutine write_step(ice,vecs,filename,time)

        implicit none 
        
        type(icesheet),         intent(IN) :: ice
        type(icesheet_vectors), intent(IN) :: vecs
        character(len=*), intent(IN) :: filename
        real(prec),       intent(IN) :: time

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
        call nc_write(filename,"T_ice",  vecs%T_ice,  units="degC",   long_name="Ice temperature",           dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"T_pmp",  vecs%T_pmp,  units="",       long_name="Ice pressure melting point",dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"T_prime",vecs%T_ice-vecs%T_pmp,units="degC",long_name="Ice temperature",     dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        
        call nc_write(filename,"cp",     vecs%cp,     units="J kg-1 K-1",   long_name="Ice heat capacity",       dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"kt",     vecs%kt,     units="J a-1 m-1 K-1",long_name="Ice thermal conductivity",dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"uz",     vecs%uz,     units="m a**-1",long_name="Ice vertical velocity",   dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"advecxy",vecs%advecxy,units="",       long_name="Ice horizontal advection",dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"Q_strn", vecs%Q_strn, units="",       long_name="Ice strain heating",      dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        
        ! Update variables (points) 
        call nc_write(filename,"Q_b",     ice%Q_b,units="",long_name="Basal heating",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"Q_geo",   ice%Q_geo,units="mW m**-2",long_name="Geothermal heat flux",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"T_srf",   ice%T_srf,units="degC",long_name="Surface temperature",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"H_ice",   ice%H_ice,units="m",long_name="Ice thickness",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"H_w",     ice%H_w,units="m",long_name="Basal water thickness",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"is_float",ice%is_float,units="",long_name="Floating flag",dim1="time",start=[n],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step

    function calc_zeta_ac(zeta_aa) result(zeta_ac)
        ! Calculate the vertical layer-edge axis (vertical ac-nodes)
        ! given the vertical layer-boundary axis as input 

        implicit none 

        real(prec), intent(IN)  :: zeta_aa(:) 
        real(prec) :: zeta_ac(size(zeta_aa)-1)

        ! Local variables
        integer :: k, nz_aa, nz_ac 

        nz_aa  = size(zeta_aa)
        nz_ac  = nz_aa - 1 

        ! Get zeta_ac (edges of zeta layers - between zeta_aa values)
        zeta_ac(1) = 0.0 
        do k = 2, nz_ac-1
            zeta_ac(k) = 0.5 * (zeta_aa(k)+zeta_aa(k+1))
        end do 
        zeta_ac(nz_ac) = 1.0 

        return 

    end function calc_zeta_ac

end program test_icetemp 
