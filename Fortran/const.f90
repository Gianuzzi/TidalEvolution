module const    
    ! Core
    integer*1, parameter :: Npar = 6 ! Amount of parameters

    ! Physical constants
    real*8, parameter  :: G      = 0.01720209895**2 ! [AU³ days⁻² Ms⁻¹]
    real*8, parameter  :: PI     = 4. * atan(1.)
    real*8, parameter  :: TWOPI  = 2. * PI  
    real*8, parameter  :: MJ2MS  = 9.54e-4         ! [Ms]
    real*8, parameter  :: MT2MS  = 3.04e-6         ! [Ms]
    real*8, parameter  :: KM2UA  = 6.68459e-9      ! [AU]
    real*8, parameter  :: YR2DAY = 365.2563        ! [day]

    ! Run
    integer*8 :: n_iter                         ! N_iterations
    integer   :: n_points                       ! N_outputs
    real*8    :: t0, tf, dt                     ! User defined constants [times]
    real*8    :: Logt, t_add                    ! Calculated time constant
    real*8    :: m0, m1                         ! Masses (constants)
    real*8    :: mu, m1p                        ! Calculated mass related constants
    real*8    :: radius0, alpha0, Q0, C0, K0cte ! Body0 constant parameters
    real*8    :: radius1, alpha1, Q1, C1, K1cte ! Body0 constant parameters

    !Extra
    character(30) :: filename    ! Output file
    character(1)  :: selection   ! In case the file exists
    logical       :: file_exists ! File existance  
    real*8        :: start_time  ! Excecution start time
    real*8        :: final_time  ! Total elapsed time
end module const