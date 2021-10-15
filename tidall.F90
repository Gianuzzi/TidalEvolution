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
    real*8 :: cos0, cos1, sin0, sin1         ! cos (o_i) & sin (o_i)
    real*8 :: fe1, fe2, fe3, fe4, fe5        ! f_i (e)

    real*8, dimension(6) :: y, ynew          ! Parameters array, P/a for iteration    
    
    contains

    ! COMMON

        real*8 function ni (a) result (nr)
            implicit none
            real*8, intent(in) :: a
            
            nr = sqrt (mu / (a**3))
        end function ni

        real*8 function Ki (a, m, r, Q) result (kr)
            implicit none
            real*8, intent(in) :: a, m, r, Q

            kr = 4.5 * G * m**2 * r**5 / (Q * ni(a))
        end function Ki
        
        real*8 function AngMom (a, e) result (AM)
            implicit none
            real*8, intent(in) :: a, e

            AM = sqrt (mu * a * (1 - e**2))
        end function AngMom
    
    ! SETTERS

        subroutine set_y (a1, e1, s1, o1, s0, o0, y)
            implicit none
            real*8, intent(in)                  :: a1, e1, s1, o1, s0, o0
            real*8, dimension(6), intent(inout) :: y
            y(1) = a1
            y(2) = e1
            y(3) = s1
            y(4) = o1
            y(5) = s0
            y(6) = o0
        end subroutine set_y

        subroutine set_fl (e, f1, f2, f3, f4, f5)
            implicit none
            real*8, intent(in)  :: e
            real*8, intent(out) :: f1, f2, f3, f4, f5
            f1 = f1_l (e)
            f2 = f2_l (e)
            f3 = f3_l (e)
            f4 = f4_l (e)
            f5 = f5_l (e)
        end subroutine set_fl

        subroutine set_fs (e, f1, f2, f3, f4, f5)
            implicit none
            real*8, intent(in)  :: e
            real*8, intent(out) :: f1, f2, f3, f4, f5
            f1 = f1_s (e)
            f2 = f2_s (e)
            f3 = f3_s (e)
            f4 = f4_s (e)
            f5 = f5_s (e)
        end subroutine set_fs

        subroutine set_cosinos (o0, o1, cos0, cos1, sin0, sin1)
            implicit none
            real*8, intent(in)  :: o0, o1
            real*8, intent(out) :: cos0, cos1, sin0, sin1
            cos0 = cos (o0)
            cos1 = cos (o1)
            sin0 = sin (o0)
            sin1 = sin (o1)
        end subroutine set_cosinos

    ! DERIVATE
        
        ! dydt (t, y__) = ynew__
        function dydtidal (t, y) result (ynew)
            implicit none
            real*8, intent(in)               :: t
            real*8, dimension(:), intent(in) :: y
            real*8, dimension(size (y))      :: ynew
            ! Here must be every dydt_i defined explicitly

            a1 = y(1)
            e1 = y(2)
            s1 = y(3)
            o1 = y(4)
            s0 = y(5)
            o0 = y(6)

            n1 = ni (a1)
            AM = AngMom(a1, e1)
            call set_fl (e1, fe1, fe2, fe3, fe4, fe5) ! Can be "set_fs" [short] if desired
            call set_cosinos (o0, o1, cos0, cos1, sin0, sin1)
            
            
            ynew(1) = 2. / m1p * a1**(-7) * & ! a1
                    & ((fe2 / n1) * (k0 * cos0 * s0 + k1 * cos1 * s1) - &
                    & fe3 * (k0 + k1))
            ynew(2) = 9. / m1p * a1**(-8) * e1 * & ! e1
                    & ((11. / 18. * fe4 / n1) * (k0 * cos0 * s0 + k1 * cos1 * s1) - &
                    & fe5 * (k0 + k1))
            ynew(3) = - k1 * n1 / (c1 * a1**6) * & !s1
                    & (fe1 * 0.5 * (1. + cos1**2) * s1 / n1 - &
                    & fe2 * cos1)
            ynew(4) = k1 * n1 / (c1 * s1 * a1**6) * & !o1
                    & (fe1 * 0.5 * s1 / n1 * (cos1 - (c1 * s1) / (m0 * AM)) -  &
                    & fe2) * sin1
            ynew(5) = - k0 * n1 / (c0 * a1**6) * & ! s0
                    & (fe1 * 0.5 * (1. + cos0**2) * s0 / n1 - &
                    & fe2 * cos0)
            ynew(6) = k0 * n1 / (c0 * s0 * a1**6) * & ! o0
                    & (fe1 * 0.5 * s0 / n1 * (cos0 - (c0 * s0) / (m1 * AM)) -  &
                    & fe2) * sin0
        end function dydtidal
end module tidall