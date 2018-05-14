program test_icetemp 

    use ncio  
    
    use defs
    use thermodynamics 
    use icetemp_grisli 
    
    implicit none 




    type icesheet 
        real(prec) :: H_ice     ! [m] Ice thickness 
        real(prec) :: H_w       ! [m] Water present at the ice base 
        real(prec) :: T_srf     ! [degC] Ice surface temperature 
        real(prec) :: Q_geo     ! [mW m-2] Geothermal heat flux 
        real(prec) :: Q_b       ! [] Basal heat production 
        logical    :: is_float  ! [-] Floating flag 
        integer    :: ibase     ! [-] Basal state (1: frozen, 2: temperate)       
        
        real(prec), allocatable :: sigma(:)  ! [-] Sigma coordinates (0: base, 1: surface)
        real(prec), allocatable :: zeta(:)   ! [-] Zeta coordinates (1:  base, 0: surface)
        
        real(prec), allocatable :: T_ice(:)   ! [degC] Ice temperature 
        real(prec), allocatable :: T_rock(:)  ! [degC] Bedrock temperature 
        real(prec), allocatable :: T_pmp(:)   ! [degC] Ice pressure melting point 
        real(prec), allocatable :: cp(:)      ! [] Ice heat capacity 
        real(prec), allocatable :: kt(:)      ! [] Ice conductivity  
        real(prec), allocatable :: uz(:)      ! [] Vertical velocity 
        real(prec), allocatable :: advecxy(:) ! [] Horizontal heat advection magnitude
        
        real(prec), allocatable :: Q_strn(:)  ! [] Strain heating 
        

    end type 

    type(icesheet) :: ice1
    
    ! Timing
    real(prec) :: t_start, t_end, dt, time  
    integer    :: n, ntot 

    t_start = 0.0     ! [yr]
    t_end   = 1000.0  ! [yr]
    dt      = 5.0     ! [yr]

    ! Calculate number of time steps to iterate and initialize time  
    ntot = 1 ! (t_end-t_start)/dt 
    time = t_start 

    ! Initialize icesheet object 
    call icesheet_allocate(ice1,nz=21,nzr=11)
    
    ! Prescribe initial eismint conditions for testing 
    call init_eismint_summit(ice1)

    do n = 1, ntot 

        ! Get current time 
        time = t_start + n*dt 

        call calc_icetemp_grisli_column(ice1%ibase,ice1%T_ice,ice1%T_rock,ice1%T_pmp, &
                                        ice1%cp,ice1%kt,ice1%uz,ice1%Q_strn,ice1%advecxy,ice1%Q_b, &
                                        ice1%Q_geo,ice1%T_srf,ice1%H_ice,ice1%H_w,ice1%is_float,dt)

    end do 

contains 


    subroutine init_eismint_summit(ice)

        implicit none 

        type(icesheet), intent(INOUT) :: ice

        ! Local variables 
        integer :: nz, nzr 

        nz  = size(ice%T_ice)
        nzr = size(ice%T_rock) 

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

        ! Make sure all vectors are deallocated
        if (allocated(ice%sigma))   deallocate(ice%sigma)
        if (allocated(ice%zeta))    deallocate(ice%zeta)
        if (allocated(ice%T_ice))   deallocate(ice%T_ice)
        if (allocated(ice%T_rock))  deallocate(ice%T_rock)
        if (allocated(ice%cp))      deallocate(ice%cp)
        if (allocated(ice%kt))      deallocate(ice%kt)
        if (allocated(ice%uz))      deallocate(ice%uz)
        if (allocated(ice%advecxy)) deallocate(ice%advecxy)
        if (allocated(ice%Q_strn))  deallocate(ice%Q_strn)
        
        ! Allocate vectors with desired lengths
        allocate(ice%sigma(nz))
        allocate(ice%zeta(nz))
        allocate(ice%T_ice(nz))
        allocate(ice%T_rock(nzr))
        allocate(ice%cp(nz))
        allocate(ice%kt(nz))
        allocate(ice%uz(nz))
        allocate(ice%advecxy(nz))
        allocate(ice%Q_strn(nz))

        ! Initialize sigma and zeta 
        ice%sigma = 0.0  
        do k = 1, nz 
            ice%sigma(k) = real(k-1,prec) / real(nz-1,prec)
        end do 
        
        ice%zeta = 1.0 - ice%sigma 

        ! Initialize remaining vectors to zero 
        ice%T_ice   = 0.0 
        ice%T_rock  = 0.0
        ice%cp      = 0.0
        ice%kt      = 0.0
        ice%uz      = 0.0
        ice%advecxy = 0.0  
        ice%Q_strn  = 0.0

        return 

    end subroutine icesheet_allocate 

    
end program test_icetemp 
