module tidall

    use const
    use f_funcs

    implicit none
    real*8 :: m0, m1, m2                     ! Masses (constants)
    real*8 :: mu, m1p, m2p                   ! Calculated mass related constants
    real*8 :: radius0, alpha0, Q0            ! Body 0 constant parameters
    real*8 :: radius1, alpha1, Q1            ! Body 1 constant parameters
    real*8 :: radius2, alpha2, Q2            ! Body 2 constant parameters
    real*8 :: C0, K01, K02                   ! Body 0 derived parameter
    real*8 :: C1, K10, K12                   ! Body 1 derived parameter
    real*8 :: C1, K20, K21                   ! Body 1 derived parameter
    real*8 :: AM, n1, n2                     ! Ang Mom, Mean Mov
    real*8 :: s0, o0                         ! Body 0 parameters
    real*8 :: a1, e1, s1, o1, w1             ! Body 1 parameters
    real*8 :: a2, e2, s2, o2  w2             ! Body 2 parameters
    real*8 :: cos0, cos1, cos2               ! cos (o_i)
    real*8 :: sin0, sin1, sin2               ! sin (o_i)
    real*8 :: fe1, fe2, fe3, fe4, fe5        ! f_i (e)

!     real*8, dimension(6) :: y, ynew          ! Parameters array, P/a for iteration
    real*8, dimension(18) :: y, ynew          ! Parameters array, P/a for iteration
    ! s0, o0, K01, K02, a1, e1, s1, o1, w1, K10, K12, 
    
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
    
        subroutine set_y0 (s0, o0, y0)
            implicit none
            real*8, intent(in)                :: s0, o0
            real*8, dimension(2), intent(out) :: y0
            
            y0(1) = s0
            y0(2) = o0
        end subroutine set_yi

        subroutine set_yi (ai, ei, si, oi, wi, yi)
            implicit none
            real*8, intent(in)                :: ai, ei, si, oi, wi
            real*8, dimension(5), intent(out) :: yi
            
            yi(1) = ai
            yi(2) = ei
            yi(3) = si
            yi(4) = oi
            yi(5) = wi
        end subroutine set_yi
        
        subroutine set_big_y (y0, y1, y2, y)
            implicit none
            real*8, dimension(2), intent(in)   :: y0
            real*8, dimension(5), intent(in)   :: y1, y2
            real*8, dimension(12), intent(out) :: y
            
            y(1: 2) = y0
            y(3: 7) = y1
            y(8:12) = y2
        end subroutine set_big_y

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

        subroutine set_cosino (oi)
            implicit none
            real*8, intent(in)  :: oi
            real*8, intent(out) :: cosi, seni
            
            cosi = cos (oi)
            sini = sin (oi)
        end subroutine set_cosino

    ! GETTERS
    
        subroutine get_from_y0 (y0, s0, o0)
            implicit none
            real*8, intent(out)              :: s0, o0
            real*8, dimension(2), intent(in) :: y0
            
            s0  = y0(1)
            o0  = y0(2)
        end subroutine get_from_y0
    
        subroutine get_from_yi (yi, ai, ei, si, oi, wi)
            implicit none
            real*8, intent(out)              :: ai, ei, si, oi, wi
            real*8, dimension(5), intent(in) :: yi
            
            ai  = yi(1) 
            ei  = yi(2) 
            si  = yi(3) 
            oi  = yi(4) 
            wi  = yi(5)
        end subroutine get_from_yi
        
        subroutine get_from_big_y (y, y0, y1, y2)
            implicit none
            real*8, dimension(12), intent(in) :: y
            real*8, dimension(2), intent(out) :: y0
            real*8, dimension(5), intent(out) :: y1, y2
            
            y0 = y(1: 2)
            y1 = y(3: 7)
            y2 = y(9:12)
        end subroutine get_from_big_y
            
        subroutine get_from_y01tidal (y01tidal, a1, e1, s1, o1, s0, o0)
            implicit none
            real*8, intent(out)              :: a1, e1, s1, o1, s0, o0
            real*8, dimension(6), intent(in) :: y01tidal
            
            a1  = y01tidal(1) 
            e1  = y01tidal(2) 
            s1  = y01tidal(3) 
            o1  = y01tidal(4) 
            s0  = y01tidal(5)
            o0  = y01tidal(6) 
        end subroutine get_from_y01tidal
    
        
    ! DERIVATE
        
        ! dydtidal0 (y0__, y1__) = ynew__
        function dydtidal_01 (m0, k0, c0, m1, k1, c1, m1p, y0, y1) result (ynew)
            implicit none
            real*8, intent(in)               :: m0, k0, c0, m1, k1, c1, m1p
            real*8, dimension(2), intent(in) :: y0
            real*8, dimension(5), intent(in) :: y1
            real*8                           :: w1 ! Dummy
            real*8, dimension(6)             :: ynew !a1, e1, s1, o1, s0, o0
            
            get_from_y0 (y0, s0, o0)
            get_from_yi (y1, a1, e1, s1, o1, w1)
            
            n1 = ni (a1)
            AM = AngMom(a1, e1)
            call set_fl (e1, fe1, fe2, fe3, fe4, fe5) ! Can be "set_fs" [short] if desired
            call set_cosino (o0, cos0, sin0)
            call set_cosino (o1, cos1, sin1)            
            
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
        end function dydtidal_01
        
        function dydtidall (t, y) result (ynew)
            implicit none
            real*8, intent(in)                 :: t
            real*8, dimension(12), intent(in)  :: y
            real*8, dimension(12), intent(out) :: ynew
            real*8, dimension(2)               :: y0
            real*8, dimension(5)               :: y1, y2
            real*8, dimension(6)               :: y01t, y02t, y12t
            
            get_from_big_y (y, y0, y1, y2)
            
            y01t = dydtidal_01 (m0, k01, c0, m1, k10, c1, m1p, y0, y1)
            y02t = dydtidal_01 (m0, k02, c0, m2, k20, c2, m2p, y0, y2)
        end function dydtidall
end module tidall
