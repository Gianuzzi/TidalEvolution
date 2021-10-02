module derivates
    use vars
    use const
    use funcs_f
    use params_f
    contains
        real*8 function dadt (a) result (der)
        implicit none
        real*8, intent(in) :: a
        real*8             :: k0a, k1a
        k0a = K0f (a)
        k1a = K1f (a)
        der = 2. / m1p * a**(-7)  * &
            & ((f2 (e1) / nf (a)) * (k0a * cos (o0) * s0 + k1a * cos (o1) * s1) - &
            & f3 (e1) * (k0a + k1a))
        end function dadt

        real*8 function dedt (e) result (der)
        implicit none
        real*8, intent(in) :: e
        der = 9. / m1p * a1**(-8) * e * &
            & ((11. / 18. * f4 (e) / n1) * (k0 * cos (o0) * s0 + k1 * cos (o1) * s1) - &
            & f5 (e) * (k0 + k1))
        end function dedt

        real*8 function ds0dt (s) result (der)
        implicit none
        real*8, intent(in) :: s
        der = - k0 * n1 / (c0 * a1**6) * &
            & (f1 (e1) * 0.5 * (1. + cos (o0)**2) * s / n1 - &
            & f2 (e1) * cos (o0))
        end function ds0dt

        real*8 function ds1dt (s) result (der)
        implicit none
        real*8, intent(in) :: s
        der = - k1 * n1 / (c1 * a1**6) * &
            & (f1 (e1) * 0.5 * (1. + cos (o1)**2) * s / n1 - &
            & f2 (e1) * cos (o1))
        end function ds1dt

        real*8 function do0dt (o) result (der)
        implicit none
        real*8, intent(in) :: o
        der = k0 * n1 / (c0 * s0 * a1**6) * &
        & (f1 (e1) * 0.5 * s0 / n1 * (cos (o) - (c0 * s0) / (m1 * AngMom (a1, e1))) -  &
        & f2 (e1)) * sin (o);
        end function do0dt

        real*8 function do1dt (o) result (der)
        implicit none
        real*8, intent(in) :: o
        der = k1 * n1 / (c1 * s1 * a1**6) * &
        & (f1 (e1) * 0.5 * s1 / n1 * (cos (o) - (c1 * s1) / (m0 * AngMom (a1, e1))) -  &
        & f2 (e1)) * sin (o);
        end function do1dt
end module derivates