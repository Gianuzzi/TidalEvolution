module tidall
    use const
    use f_funcs

    implicit none
    real*8 :: m0, m1                         ! Masses (constants)
    real*8 :: mu, m1p                        ! Calculated mass related constants
    real*8 :: radius0, alpha0, Q0            ! Body 0 constant parameters
    real*8 :: radius1, alpha1, Q1            ! Body 1 constant parameters
    real*8 :: C0, K0                         ! Body 0 derived parameter
    real*8 :: C1, K1                         ! Body 1 derived parameter
    real*8 :: AM, n1                         ! Ang Mom, Mean Mov
    real*8 :: s0, o0                         ! Body 0 parameters
    real*8 :: a1, e1, s1, o1                 ! Body 1 parameters
    real*8 :: a10, e10, s10, o10, o00, s00   ! Iteration variables
    real*8 :: cos0, cos1                     ! cos (o_i)
    real*8 :: fe1, fe2, fe3, fe4, fe5        ! f_i (e)
    
    contains

    ! COMMON

        real*8 function ni (a) result (nr)
            implicit none
            real*8, intent (in) :: a
            
            nr = sqrt (mu / (a**3))
        end function ni

        real*8 function Ki (a, m, r, Q) result (kr)
            implicit none
            real*8, intent (in) :: a, m, r, Q

            kr = 4.5 * G * m**2 * r**5 / (Q * ni(a))
        end function Ki
        
        real*8 function AngMom (a, e) result (AM)
            implicit none
            real*8, intent (in) :: a, e

            AM = sqrt (mu * a * (1 - e**2))
        end function AngMom
    
    ! SETTERS

        subroutine set_fs (e, f1r, f2r, f3r, f4r, f5r)
            implicit none
            real*8, intent (in)  :: e
            real*8, intent (out) :: f1r, f2r, f3r, f4r, f5r
            f1r = f1_l (e)
            f2r = f2_l (e)
            f3r = f3_l (e)
            f4r = f4_l (e)
            f5r = f5_l (e)
        end subroutine set_fs

        subroutine set_fshorts (e, f1r, f2r, f3r, f4r, f5r)
            implicit none
            real*8, intent (in)  :: e
            real*8, intent (out) :: f1r, f2r, f3r, f4r, f5r
            f1r = f1_s (e)
            f2r = f2_s (e)
            f3r = f3_s (e)
            f4r = f4_s (e)
            f5r = f5_s (e)
        end subroutine set_fshorts

        subroutine set_cosos (o0i, o1i, coso0, coso1)
            implicit none
            real*8, intent (in)  :: o0i, o1i
            real*8, intent (out) :: coso0, coso1
            coso0 = cos (o0i)
            coso1 = cos (o1i)
        end subroutine set_cosos

    ! DERIVATES

        real*8 function dadt (t, a) result (der)
            implicit none
            real*8, intent(in) :: t, a

            der = 2. / m1p * a**(-7)  * &
                & ((fe2 / ni (a)) * (k0 * cos0 * s0 + k1 * cos1 * s1) - &
                & fe3 * (k0 + k1))
        end function dadt

        real*8 function dedt (t, e) result (der)
            implicit none
            real*8, intent(in) :: t, e

            der = 9. / m1p * a1**(-8) * e * &
                & ((11. / 18. * f4_l (e) / n1) * (k0 * cos0 * s0 + k1 * cos1 * s1) - &
                & f5_l (e) * (k0 + k1))
        end function dedt

        real*8 function ds0dt (t, s) result (der)
            implicit none
            real*8, intent(in) :: t, s

            der = - k0 * n1 / (c0 * a1**6) * &
                & (fe1 * 0.5 * (1. + cos0**2) * s / n1 - &
                & fe2 * cos0)
        end function ds0dt

        real*8 function ds1dt (t, s) result (der)
            implicit none
            real*8, intent(in) :: t, s

            der = - k1 * n1 / (c1 * a1**6) * &
                & (fe1 * 0.5 * (1. + cos1**2) * s / n1 - &
                & fe2 * cos1)
        end function ds1dt

        real*8 function do0dt (t, o) result (der)
            implicit none
            real*8, intent(in) :: t, o

            der = k0 * n1 / (c0 * s0 * a1**6) * &
            & (fe1 * 0.5 * s0 / n1 * (cos (o) - (c0 * s0) / (m1 * AM)) -  &
            & fe2) * sin (o);
        end function do0dt

        real*8 function do1dt (t, o) result (der)
            implicit none
            real*8, intent(in) :: t, o

            der = k1 * n1 / (c1 * s1 * a1**6) * &
            & (fe1 * 0.5 * s1 / n1 * (cos (o) - (c1 * s1) / (m0 * AM)) -  &
            & fe2) * sin (o);
        end function do1dt

        real*8 function dedt_s (t, e) result (der)
            implicit none
            real*8, intent(in) :: t, e

            der = 9. / m1p * a1**(-8) * e * &
                & ((11. / 18. * f4_s (e) / n1) * (k0 * cos0 * s0 + k1 * cos1 * s1) - &
                & f5_s (e) * (k0 + k1))
        end function dedt_s

end module tidall