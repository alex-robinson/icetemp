program test_icetemp 

    use ncio  
    
    use defs
    use thermodynamics 
    use icetemp_grisli 
    
    implicit none 



    type icesheet_vectors
        real(prec), allocatable :: sigma(:)   ! [-] Sigma coordinates (either height or depth)
        real(prec), allocatable :: T_ice(:)   ! [degC] Ice temperature 
        real(prec), allocatable :: T_rock(:)  ! [degC] Bedrock temperature 
        real(prec), allocatable :: T_pmp(:)   ! [degC] Ice pressure melting point 
        real(prec), allocatable :: cp(:)      ! [] Ice heat capacity 
        real(prec), allocatable :: kt(:)      ! [] Ice conductivity  
        real(prec), allocatable :: uz(:)      ! [] Vertical velocity 
        real(prec), allocatable :: advecxy(:) ! [] Horizontal heat advection magnitude
        real(prec), allocatable :: Q_strn(:)  ! [] Strain heating 
        
    end type 

    type icesheet 
        real(prec) :: H_ice     ! [m] Ice thickness 
        real(prec) :: H_w       ! [m] Water present at the ice base 
        real(prec) :: T_srf     ! [degC] Ice surface temperature 
        real(prec) :: smb       ! [m a**-1] Surface mass balance
        real(prec) :: bmb       ! [m a**-1] Basal mass balance
        real(prec) :: Q_geo     ! [mW m-2] Geothermal heat flux 
        real(prec) :: Q_b       ! [] Basal heat production 
        logical    :: is_float  ! [-] Floating flag 
        integer    :: ibase     ! [-] Basal state (1: frozen, 2: temperate)       
        
        type(icesheet_vectors) :: up     ! For height coordinate systems with k=1 base and k=nz surface
        type(icesheet_vectors) :: dwn    ! For depth coordinate systems with k=1 surface and k=nz base 
        
    end type 

    type(icesheet) :: ice1
    
    ! Timing
    real(prec) :: t_start, t_end, dt, time  
    integer    :: n, ntot 

    ! Output 
    character(len=512) :: file1D 
    real(prec)         :: dt_out 

    t_start = 0.0     ! [yr]
    t_end   = 10000.0  ! [yr]
    dt      = 5.0     ! [yr]

    file1D  = "test.nc" 
    dt_out  = 100.0      ! [yr] 

    ! Calculate number of time steps to iterate and initialize time  
    ntot = (t_end-t_start)/dt 
    time = t_start 

    ! Initialize icesheet object 
    call icesheet_allocate(ice1,nz=21,nzr=11)
    
    ! Prescribe initial eismint conditions for testing 
    call init_eismint_summit(ice1)

    ! Initialize output file and write intial conditions 
    call write_init(ice1,filename=file1D,sigma=ice1%up%sigma,time_init=time)
    call write_step(ice1,ice1%up,filename=file1D,time=time)

    ! Transfer info to ice1%dwn 
    call icesheet_up_to_dwn(ice1%dwn,ice1%up)

    do n = 1, ntot 

        ! Get current time 
        time = t_start + n*dt 

        ! Call grisli icetemp routine
        call calc_icetemp_grisli_column(ice1%ibase,ice1%dwn%T_ice,ice1%dwn%T_rock,ice1%dwn%T_pmp, &
                                        ice1%dwn%cp,ice1%dwn%kt,ice1%dwn%uz,ice1%dwn%Q_strn,ice1%dwn%advecxy, &
                                        ice1%Q_b,ice1%Q_geo,ice1%T_srf,ice1%H_ice,ice1%H_w,ice1%is_float,dt)
    
        ! Transfer info back to ice1%up for writing
        call icesheet_dwn_to_up(ice1%up,ice1%dwn)

        if (mod(time,dt_out)==0) then 
            call write_step(ice1,ice1%up,filename=file1D,time=time)
        end if 

        if (mod(time,20.0)==0) then
            write(*,"(a,f14.4)") "time = ", time
        end if 

    end do 

contains 


    subroutine init_eismint_summit(ice)

        implicit none 

        type(icesheet), intent(INOUT) :: ice

        ! Local variables 
        integer :: k, nz, nzr 

        nz  = size(ice%up%T_ice)
        nzr = size(ice%up%T_rock) 

        ice%T_srf    = 239.0-T0    ! [degC]
        ice%smb      = 0.5         ! [m/a]
        ice%bmb      = 0.0         ! [m/a]
        ice%Q_geo    = 42.0        ! [mW/m2]
        ice%H_ice    = 2997.0      ! [m] Summit thickness
        ice%H_w      = 0.0         ! [m] No basal water
        ice%Q_b      = 0.0         ! [] No basal frictional heating 
        ice%is_float = .FALSE.     ! Grounded point 
        ice%ibase    = 1           ! Frozen 

        ice%up%Q_strn  = 0.0         ! [] No internal strain heating 
        ice%up%advecxy = 0.0         ! [] No horizontal advection 

        ! Calculate pressure melting point 
        ice%up%T_pmp = calc_T_pmp((1.0-ice%up%sigma)*ice%H_ice,T0) - T0 

        ! Define surface temperature of ice based on simple atmospheric correction below zero
        ice%up%T_ice(nz) = ice%T_srf 

        ! Basal temperature is 10 deg below freezing point 
        ice%up%T_ice(1) = ice%up%T_pmp(1) - 10.0 

        ! Intermediate layers are linearly interpolated 
        do k = 2, nz-1 
            ice%up%T_ice(k) = ice%up%T_ice(1)+ice%up%sigma(k)*(ice%up%T_ice(nz)-ice%up%T_ice(1))
        end do 

        ice%up%T_rock = ice%up%T_ice(1) 

        ! Define vertical velocity profile (linear)
        ice%up%uz(nz) = -ice%smb 
        ice%up%uz(1)  = 0.0 
        do k = 2, nz-1 
            ice%up%uz(k) = ice%up%uz(1)+ice%up%sigma(k)*(ice%up%uz(nz)-ice%up%uz(1))
        end do 

        return 

    end subroutine init_eismint_summit 

    subroutine icesheet_allocate(ice,nz,nzr)
        ! Allocate the ice sheet object 

        implicit none 

        type(icesheet), intent(INOUT) :: ice 
        integer, intent(IN) :: nz            ! Number of ice points
        integer, intent(IN) :: nzr           ! Number of rock points 

        ! Local variables 
        integer :: k 

        ! First allocate 'up' variables (with vertical coordinate as height)

        ! Make sure all vectors are deallocated
        if (allocated(ice%up%sigma))   deallocate(ice%up%sigma)
        if (allocated(ice%up%T_ice))   deallocate(ice%up%T_ice)
        if (allocated(ice%up%T_rock))  deallocate(ice%up%T_rock)
        if (allocated(ice%up%T_pmp))   deallocate(ice%up%T_pmp)
        if (allocated(ice%up%cp))      deallocate(ice%up%cp)
        if (allocated(ice%up%kt))      deallocate(ice%up%kt)
        if (allocated(ice%up%uz))      deallocate(ice%up%uz)
        if (allocated(ice%up%advecxy)) deallocate(ice%up%advecxy)
        if (allocated(ice%up%Q_strn))  deallocate(ice%up%Q_strn)
        
        ! Allocate vectors with desired lengths
        allocate(ice%up%sigma(nz))
        allocate(ice%up%T_ice(nz))
        allocate(ice%up%T_rock(nzr))
        allocate(ice%up%T_pmp(nz))
        allocate(ice%up%cp(nz))
        allocate(ice%up%kt(nz))
        allocate(ice%up%uz(nz))
        allocate(ice%up%advecxy(nz))
        allocate(ice%up%Q_strn(nz))

        ! Initialize sigma and zeta 
        ice%up%sigma = 0.0  
        do k = 1, nz 
            ice%up%sigma(k) = real(k-1,prec) / real(nz-1,prec)
        end do 

        ! Initialize remaining vectors to zero 
        ice%up%T_ice   = 0.0 
        ice%up%T_rock  = 0.0
        ice%up%T_pmp   = 0.0 
        ice%up%cp      = 0.0
        ice%up%kt      = 0.0
        ice%up%uz      = 0.0
        ice%up%advecxy = 0.0  
        ice%up%Q_strn  = 0.0

        ! Now allocate down variables too (with vertical coordinate as depth)
        ice%dwn = ice%up 
        ice%dwn%sigma = 1.0 - ice%up%sigma 

        write(*,*) "Allocated icesheet variables."

        return 

    end subroutine icesheet_allocate 

    subroutine write_init(ice,filename,sigma,time_init)

        implicit none 

        type(icesheet),   intent(IN) :: ice 
        character(len=*), intent(IN) :: filename 
        real(prec),       intent(IN) :: sigma(:) 
        real(prec),       intent(IN) :: time_init

        ! Local variables 
        integer    :: nzr 
        real(prec) :: dzr 

        nzr = size(ice%up%T_rock,1)

        ! Initialize netcdf file and dimensions
        call nc_create(filename)
        call nc_write_dim(filename,"sigma", x=sigma, units="1")
        call nc_write_dim(filename,"sigmar",x=0,nx=nzr,dx=1,units="1")
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
        character(len=12), parameter :: vert_dim = "sigma"

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
        if (abs(time-time_prev).gt.1e-5) n = n+1 

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! Update variables (vectors) 
        call nc_write(filename,"T_ice",  vecs%T_ice,  units="degC",   long_name="Ice temperature",         dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"T_rock", vecs%T_rock, units="degC",   long_name="Bedrock temperature",     dim1="sigmar",dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"T_pmp",  vecs%T_pmp,  units="",       long_name="Ice pressure melting point",dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"cp",     vecs%cp,     units="",       long_name="Ice heat capacity",       dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"kt",     vecs%kt,     units="",       long_name="Ice thermal conductivity",dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"uz",     vecs%uz,     units="m a**-1",long_name="Ice vertical velocity",   dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"advecxy",vecs%advecxy,units="",       long_name="Ice horizontal advection",dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"Q_strn", vecs%Q_strn, units="",       long_name="Ice strain heating",      dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        
        ! Update variables (points) 
        call nc_write(filename,"ibase",   ice%ibase,units="",long_name="Basal state",dim1="time",start=[n],ncid=ncid)
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

    subroutine icesheet_up_to_dwn(dwn,up)

        implicit none 

        type(icesheet_vectors), intent(INOUT) :: dwn 
        type(icesheet_vectors), intent(IN)    :: up 

        ! Local variables 
        integer :: k, nz, nzr 

        nz  = size(up%T_ice,1)
        nzr = size(up%T_rock,1)

        do k = 1, nz 
            dwn%T_ice(nz-k+1)   = up%T_ice(k) 
            dwn%T_pmp(nz-k+1)   = up%T_pmp(k) 
            dwn%cp(nz-k+1)      = up%cp(k) 
            dwn%kt(nz-k+1)      = up%kt(k) 
            dwn%uz(nz-k+1)      = up%uz(k) 
            dwn%advecxy(nz-k+1) = up%advecxy(k) 
            dwn%Q_strn(nz-k+1)  = up%Q_strn(k) 
        end do 

        dwn%T_rock = up%T_rock 

        return 

    end subroutine icesheet_up_to_dwn 

    subroutine icesheet_dwn_to_up(up,dwn)

        implicit none 

        type(icesheet_vectors), intent(INOUT) :: up 
        type(icesheet_vectors), intent(IN)    :: dwn 

        ! Local variables 
        integer :: k, nz, nzr 

        nz  = size(dwn%T_ice,1)
        nzr = size(dwn%T_rock,1)

        do k = 1, nz 
            up%T_ice(nz-k+1)   = dwn%T_ice(k) 
            up%T_pmp(nz-k+1)   = dwn%T_pmp(k) 
            up%cp(nz-k+1)      = dwn%cp(k) 
            up%kt(nz-k+1)      = dwn%kt(k) 
            up%uz(nz-k+1)      = dwn%uz(k) 
            up%advecxy(nz-k+1) = dwn%advecxy(k) 
            up%Q_strn(nz-k+1)  = dwn%Q_strn(k) 
        end do 

        up%T_rock = dwn%T_rock 
        
        return 

    end subroutine icesheet_dwn_to_up 

end program test_icetemp 
