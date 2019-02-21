module defs
    
    implicit none 

    ! =========================================================================
    !
    ! CONSTANTS (program precision, global constants)
    !
    ! =========================================================================

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: prec = sp 

    ! Write flags 
    logical, parameter :: write_log = .TRUE. 

    ! Missing value and aliases
    real(prec), parameter :: MISSING_VALUE_DEFAULT = real(-9999.0,prec)
    real(prec), parameter :: MISSING_VALUE = MISSING_VALUE_DEFAULT
    real(prec), parameter :: MV = MISSING_VALUE_DEFAULT
    
    ! Error distance (very large), error index, and smallest number epsilon 
    real(prec), parameter :: ERR_DIST = real(1e8,prec) 
    integer,    parameter :: ERR_IND  = -1 
    real(prec), parameter :: eps      = real(1e-8,prec) 
    
    ! Mathematical constants
    real(prec), parameter :: pi  = real(2._dp*acos(0.0_dp),prec)
    real(prec), parameter :: degrees_to_radians = real(pi / 180._dp,prec)  ! Conversion factor between radians and degrees
    real(prec), parameter :: radians_to_degrees = real(180._dp / pi,prec)  ! Conversion factor between degrees and radians
    
    ! The constants below should be loaded using the global subroutine
    ! defined below `global_init`.
    ! Note: The key limitation imposed by defining the parameters defined 
    ! globally is that these constants must be the same for all domains 
    ! being run in the same program. 
    
    ! Physical constants 
    real(prec), parameter :: sec_year = 31556926.0   ! [s] EISMINT value
    real(prec), parameter :: g        = 9.81         ! [m s-2] Gravitational accel.
    real(prec), parameter :: T0       = 273.15       ! [K] Reference freezing temperature 
    real(prec), parameter :: rho_ice  =  910.0       ! [kg m-3] Density ice
    real(prec), parameter :: rho_w    = 1000.0       ! [kg m-3] Density water      
    real(prec), parameter :: rho_sw   = 1028.0       ! [kg m-3] Density seawater      
    real(prec), parameter :: rho_a    = 3300.0       ! [kg m-3] Density asthenosphere
    real(prec), parameter :: rho_m    = 2000.0       ! [kg m-3] Density mantle (lithosphere)
    real(prec), parameter :: L_ice    = 333500.0     ! [J kg-1] Latent heat

    public   ! All defs are public

contains 

end module defs

