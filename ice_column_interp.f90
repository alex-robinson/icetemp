module ice_column_interp 

    use defs, only : prec, pi
    
    implicit none 


    type ice_column_hires_type
        integer :: nz_aa, nz_ac  
        real(prec), allocatable :: enth(:)        ! nz_aa [J kg] Ice column enthalpy
        real(prec), allocatable :: T_ice(:)       ! nz_aa [K] Ice column temperature
        real(prec), allocatable :: omega(:)       ! nz_aa [-] Ice column water content fraction
        real(prec), allocatable :: T_pmp(:)       ! nz_aa [K] Pressure melting point temp.
        real(prec), allocatable :: cp(:)          ! nz_aa [J kg-1 K-1] Specific heat capacity
        real(prec), allocatable :: kt(:)          ! nz_aa [J a-1 m-1 K-1] Heat conductivity 
        real(prec), allocatable :: advecxy(:)     ! nz_aa [K a-1] Horizontal heat advection 
        real(prec), allocatable :: uz(:)          ! nz_ac [m a-1] Vertical velocity 
        real(prec), allocatable :: Q_strn(:)      ! nz_aa [J a-1 m-3] Internal strain heat production in ice
        real(prec), allocatable :: zeta_aa(:)     ! nz_aa [--] Vertical sigma coordinates (zeta==height), layer centered aa-nodes
        real(prec), allocatable :: zeta_ac(:)     ! nz_ac [--] Vertical height axis temperature (0:1), layer edges ac-nodes
        real(prec), allocatable :: dzeta_a(:)     ! nz_aa [--] Solver discretization helper variable ak
        real(prec), allocatable :: dzeta_b(:)     ! nz_aa [--] Solver discretization helper variable bk
    end type 


contains 



    subroutine ice_column_hires_init(dat,nz_aa,nz_ac,fac,zeta_scale,zeta_exp)

        implicit none 

        type(ice_column_hires_type), intent(INOUT) :: dat 
        integer,          intent(IN) :: nz_aa 
        integer,          intent(IN) :: nz_ac 
        integer,          intent(IN) :: fac 
        character(len=*), intent(IN) :: zeta_scale 
        real(prec),       intent(IN) :: zeta_exp 

        ! Local variables 
        integer :: nz_aa_hi
        integer :: nz_ac_hi 

        if (fac .lt. 1) then 
            write(*,*) "ice_column_hires_init:: Error: factor scaling only works for fac > 1."
            stop 
        end if 

        ! Calculate high-res number of grid points based on integer scaling
        nz_aa_hi = (nz_aa-1)*fac + 1
        nz_ac_hi = nz_aa - 1 

        ! Allocate the column vectors 
        call ice_column_hires_alloc(dat,nz_aa_hi,nz_ac_hi)

        ! Calculate the axis variables 
        call calc_zeta(dat%zeta_aa,dat%zeta_ac,zeta_scale,zeta_exp)
        call calc_dzeta_terms(dat%dzeta_a,dat%dzeta_b,dat%zeta_aa,dat%zeta_ac)

        return 

    end subroutine ice_column_hires_init 

    subroutine calc_zeta(zeta_aa,zeta_ac,zeta_scale,zeta_exp)
        ! Calculate the vertical layer-edge axis (vertical ac-nodes)
        ! and the vertical cell-center axis (vertical aa-nodes),
        ! including an extra zero-thickness aa-node at the base and surface

        implicit none 

        real(prec), intent(INOUT)  :: zeta_aa(:) 
        real(prec), intent(INOUT)  :: zeta_ac(:) 
        character(*), intent(IN)   :: zeta_scale 
        real(prec),   intent(IN)   :: zeta_exp 

        ! Local variables
        integer :: k, nz_aa, nz_ac 

        nz_aa  = size(zeta_aa)
        nz_ac  = size(zeta_ac)   ! == nz_aa - 1 

        ! Initially define a linear zeta scale 
        ! Base = 0.0, Surface = 1.0 
        do k = 1, nz_ac
            zeta_ac(k) = 0.0 + 1.0*(k-1)/real(nz_ac-1)
        end do 

        ! Scale zeta to produce different resolution through column if desired
        ! zeta_scale = ["linear","exp","wave"]
        select case(trim(zeta_scale))
            
            case("exp")
                ! Increase resolution at the base 
                zeta_ac = zeta_ac**(zeta_exp) 

            case("tanh")
                ! Increase resolution at base and surface 

                zeta_ac = tanh(1.0*pi*(zeta_ac-0.5))
                zeta_ac = zeta_ac - minval(zeta_ac)
                zeta_ac = zeta_ac / maxval(zeta_ac)

            case DEFAULT
            ! Do nothing, scale should be linear as defined above
        
        end select  
        
        ! Get zeta_aa (between zeta_ac values, as well as at the base and surface)
        zeta_aa(1) = 0.0 
        do k = 2, nz_aa-1
            zeta_aa(k) = 0.5 * (zeta_ac(k-1)+zeta_ac(k))
        end do 
        zeta_aa(nz_aa) = 1.0 

        return 

    end subroutine calc_zeta
    
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
    
    subroutine ice_column_hires_alloc(dat,nz_aa,nz_ac)

        implicit none 

        type(ice_column_hires_type), intent(INOUT) :: dat 
        integer, intent(IN) :: nz_aa 
        integer, intent(IN) :: nz_ac 

        call ice_column_hires_dealloc(dat)

        allocate(dat%enth(nz_aa))
        allocate(dat%T_ice(nz_aa))
        allocate(dat%T_pmp(nz_aa))
        allocate(dat%cp(nz_aa))
        allocate(dat%kt(nz_aa))
        allocate(dat%uz(nz_ac))         ! ac-nodes 
        allocate(dat%Q_strn(nz_aa))
        allocate(dat%zeta_aa(nz_aa))
        allocate(dat%zeta_ac(nz_aa))
        allocate(dat%dzeta_a(nz_aa))
        allocate(dat%dzeta_b(nz_aa))
        
        ! Store high-res npts 
        dat%nz_aa = nz_aa 
        dat%nz_ac = nz_ac 

        return 

    end subroutine ice_column_hires_alloc

    subroutine ice_column_hires_dealloc(dat)

        implicit none 

        type(ice_column_hires_type), intent(INOUT) :: dat 
        
        if (allocated(dat%enth))    deallocate(dat%enth)
        if (allocated(dat%T_ice))   deallocate(dat%T_ice)
        if (allocated(dat%T_pmp))   deallocate(dat%T_pmp)
        if (allocated(dat%cp))      deallocate(dat%cp)
        if (allocated(dat%kt))      deallocate(dat%kt)
        if (allocated(dat%advecxy)) deallocate(dat%advecxy)
        if (allocated(dat%uz))      deallocate(dat%uz)
        if (allocated(dat%Q_strn))  deallocate(dat%Q_strn)
        if (allocated(dat%zeta_aa)) deallocate(dat%zeta_aa)
        if (allocated(dat%zeta_ac)) deallocate(dat%zeta_ac)
        if (allocated(dat%dzeta_a)) deallocate(dat%dzeta_a)
        if (allocated(dat%dzeta_b)) deallocate(dat%dzeta_b)

        return 

    end subroutine ice_column_hires_dealloc
    
end module ice_column_interp 