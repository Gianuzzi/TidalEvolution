module tidal_vars
    real*8 :: s0, o0                         ! Body 0 parameters
    real*8 :: a1, e1, s1, o1                 ! Body 1 parameters
    real*8 :: m0, m1                         ! Masses (constants)
    real*8 :: mu, m1p                        ! Calculated mass related constants
    real*8 :: radius0, alpha0, Q0, C0, K0cte ! Body 0 constant parameters
    real*8 :: radius1, alpha1, Q1, C1, K1cte ! Body 1 constant parameters
    real*8 :: K0                             ! Body 0 derived parameter
    real*8 :: K1, n1                         ! Body 1 derived parameter
    real*8 :: a10, e10, s10, o10, o00, s00   ! Iteration variables

    contains
    real*8 function nf (a) result (ni)
        implicit none
        real*8, intent (in) :: a
        
        ni = sqrt (mu / (a**3))
    end function nf

    real*8 function K0f (a) result (k0r)
        implicit none
        real*8, intent (in) :: a

        k0r = K0cte * a**(1.5)
    end function K0f

    real*8 function K1f (a) result (k1r)
        implicit none
        real*8, intent (in) :: a

        k1r = K1cte * a**(1.5)
    end function K1f

    real*8 function AngMom (a, e) result (AM)
        implicit none
        real*8, intent (in) :: a, e

        AM = sqrt (mu * a * (1 - e**2))
    end function AngMom
end module tidal_vars