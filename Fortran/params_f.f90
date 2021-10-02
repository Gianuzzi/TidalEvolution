module params_f
    use const
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
end module params_f