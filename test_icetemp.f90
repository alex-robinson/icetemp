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

        real(prec), allocatable :: sigma(:)  ! [-] Sigma coordinates (0: base, 1: surface)
        real(prec), allocatable :: zeta(:)   ! [-] Zeta coordinates (1:  base, 0: surface)
        
        real(prec), allocatable :: T_ice(:)   ! [degC] Ice temperature 
        real(prec), allocatable :: T_rock(:)  ! [degC] Bedrock temperature 
        real(prec), allocatable :: Q_strn(:)  ! [] Strain heating 
        
    end type 

    type(icesheet) :: ice1

    call icesheet_allocate(ice1,nz=21,nzr=11)
    

contains 


    subroutine init_eismint_summit()

        implicit none 



        return 

    end subroutine init_eismint_summit 

    subroutine icesheet_allocate(ice,nz,nzr)
        ! Allocate the ice sheet object 

        implicit none 

        type(icesheet), intent(INOUT) :: ice 
        integer, intent(IN) :: nz            ! Number of ice points
        integer, intent(IN) :: nzr           ! Number of rock points 

        if (allocated(ice%sigma))  deallocate(ice%sigma)
        if (allocated(ice%zeta))   deallocate(ice%zeta)
        if (allocated(ice%T_ice))  deallocate(ice%T_ice)
        if (allocated(ice%T_rock)) deallocate(ice%T_rock)
        if (allocated(ice%Q_strn)) deallocate(ice%Q_strn)
        
        allocate(ice%sigma(nz))
        allocate(ice%zeta(nz))
        allocate(ice%T_ice(nz))
        allocate(ice%T_rock(nzr))
        allocate(ice%Q_strn(nz))
        

        return 

    end subroutine icesheet_allocate 

    
end program test_icetemp 
