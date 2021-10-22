module const

    implicit none    
    real*8, parameter  :: G       = 0.01720209895**2 ! [AU³ days⁻² Ms⁻¹]
    real*8, parameter  :: PI      = 4. * atan (1.)   ! [rad]
    real*8, parameter  :: TWOPI   = 2. * PI          ! [rad]
    real*8, parameter  :: DEG2RAD = PI / 180.        ! [rad]
    real*8, parameter  :: MJ2MS   = 9.54e-4          ! [Ms]
    real*8, parameter  :: MT2MS   = 3.04e-6          ! [Ms]
    real*8, parameter  :: KM2UA   = 6.68459e-9       ! [AU]
    real*8, parameter  :: RS2UA   = 4.65047e-3       ! [AU]
    real*8, parameter  :: RJ2UA   = 4.77895e-4       ! [AU]
    real*8, parameter  :: RT2UA   = 4.26352e-5       ! [AU]
    real*8, parameter  :: YR2DAY  = 365.2563         ! [day]

end module const