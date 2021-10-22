module tidall_ip

    use const

    implicit none
    integer, parameter :: N_bodies = 3

    real*8, dimension(N_bodies) :: m            ! Masses (constants)
    real*8, dimension(N_bodies) :: r, alp, C    ! Bodies constant parameters
    real*8, dimension(N_bodies) :: Q, k2dt      ! Bodies derived parameters
    real*8, dimension(N_bodies) :: mu, mp, n    ! Bodies 1... derived parameters
    
    real*8, dimension(N_bodies, N_bodies) :: TSk ! TIDAL STRENGTH [k]

    real*8 :: s0, o0                         ! Body 0 parameters
    real*8 :: a1, e1, s1, o1, vp1            ! Body 1 parameters
    real*8 :: a2, e2, s2, o2, vp2            ! Body 2 parameters

    real*8, dimension(2)  :: y0              ! Body 0 Parameters array
    real*8, dimension(5)  :: y1, y2          ! Body 1 and Body 2 Parameters array

    real*8, dimension(2 + 5 * (N_bodies - 1)) :: y, ynew  ! Parameters array, and same for iteration
    
    contains

    ! COMMON

        real*8 function mu_f (m0, mj) result (mu)
            implicit none
            real*8, intent(in) :: m0, mj
            
            mu = G * (m0 + mj)
        end function mu_f

        real*8 function mprime_f (mi, mj) result (mp)
            implicit none
            real*8, intent(in) :: mi, mj

            mp = (mi * mj) / (mi + mj)
        end function mprime_f

        real*8 function n_f (a, mu) result (n)
            implicit none
            real*8, intent(in) :: a, mu
            
            n = sqrt (mu / (a**3))
        end function n_f

        real*8 function TidalStrength_f (n, mj, radii, Q) result (k)
            implicit none
            real*8, intent(in) :: n, mj, radii, Q

            k = 4.5 * G * mj**2 * radii**5 / (Q * n)
        end function TidalStrength_f
        
        real*8 function AngMom_f (a, e, mu) result (AM)
            implicit none
            real*8, intent(in) :: a, e, mu

            AM = sqrt (mu * a * (1 - e**2))
        end function AngMom_f

        real*8 function K_f (e, vp) result (K)
            implicit none
            real*8, intent(in) :: e, vp

            K = e * cos (vp)
        end function K_f

        real*8 function H_f (e, vp) result (H)
            implicit none
            real*8, intent(in) :: e, vp

            H = e * sin (vp)
        end function H_f
    
    ! SETTERS
        subroutine set_initial_params (a, Q, alpha, radius, m, C, mu, mp, n, k2dt, K)
            implicit none
            real*8, dimension(:), intent(in)                  :: a, m, Q, alpha, radius
            real*8, dimension(size (Q)), intent(out)          :: mu, mp, n, C, k2dt
            real*8, dimension(size (Q), size(Q)), intent(out) :: K
            integer                                           :: i, j
            
            K = 0.
            do i = size (alpha), 1, -1 !Backwards, for n(1) setting
                if (i > 1) then
                    mp(i) = mprime_f (sum (m(:(i - 1))), m(i))
                    mu(i) = mu_f (sum (m(:(i - 1))), m(i))                    
                    n(i)  = n_f (a(i - 1), mu(i))
                else
                    mu(i) = 0.       ! Unused
                    mp(i) = 0.       ! Unused
                    n(i)  = n(i + 1) ! n(1) = n(2)
                    
                end if
                C(i)    = alpha(i) * m(i) * radius(i)**2
                k2dt(i) = 1.5 / Q(i) / n(i)                
                do j = 1, N_bodies
                    if (i .ne. j) then
                        K(i, j) = TidalStrength_f(n(i), m(j), radius(i), Q(i))
                    end if
                end do
            end do

        end subroutine set_initial_params

        subroutine set_cosino (o, cosi, sini)
            implicit none
            real*8, intent(in)  :: o
            real*8, intent(out) :: cosi, sini
            
            cosi = cos (o)
            sini = sin (o)
        end subroutine set_cosino

        subroutine set_y0 (s0, o0, y0)
            implicit none
            real*8, intent(in)                :: s0, o0
            real*8, dimension(2), intent(out) :: y0
            
            y0(1) = s0
            y0(2) = o0
        end subroutine set_y0

        subroutine set_yi (a, e, s, o, vp, yi)
            implicit none
            real*8, intent(in)                :: a, e, s, o, vp
            real*8, dimension(5), intent(out) :: yi
            
            yi(1) = a
            yi(2) = e
            yi(3) = s
            yi(4) = o
            yi(5) = vp
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

        subroutine set_fl(e, f1, f2, f3 ,f4, f5)
            real*8, intent(in)  :: e
            real*8, intent(out) :: f1, f2, f3, f4, f5
            real*8              :: e2, e4, e6, e8, de
            e2 = e * e
            e4 = e2 * e2
            e6 = e2 * e4
            e8 = e2 * e6
            de = (1. - e2)
            f1 = (1. + 3.0 * e2 + 0.375 * e4) / de**(4.5)
            f2 = (1. + 7.5 * e2 + 5.625 * e4 + 0.3125 * e6) / de**6
            f3 = (1. + 15.5 * e2 + 31.875 * e4 + 11.5625 * e6 + 0.390625 * e8) / de**(7.5)
            f4 = (1. + 1.5 * e2 + 0.125 * e4) / de**5
            f5 = (1. + 3.75 * e2 + 1.875 * e4 + 0.078125 * e6) / de**(6.5)
        end subroutine set_fl

        subroutine set_fs(e, f1, f2, f3 ,f4, f5)
            real*8, intent(in)  :: e
            real*8, intent(out) :: f1, f2, f3, f4, f5
            real*8              :: e2
            e2 = e * e
            f1 =  1. + 7.5 * e2
            f2 =  1. + 13.5 * e2
            f3 =  1. + 23. * e2
            f4 =  1. + 6.5 * e2
            f5 =  1. + 10.25 * e2
        end subroutine set_fs

    ! GETTERS

        subroutine get_evp(K, H, e, vp)
            implicit none
            real*8, intent(in)  :: K, H
            real*8, intent(out) :: e, vp

            e  = sqrt (K**2 + H**2)
            vp = atan2 (H, K)        
        end subroutine get_evp
            
        subroutine get_from_y0 (y0, s0, o0)
            implicit none
            real*8, intent(out)              :: s0, o0
            real*8, dimension(2), intent(in) :: y0
            
            s0 = y0(1)
            o0 = y0(2)
        end subroutine get_from_y0
    
        subroutine get_from_yi (yi, a, e, s, o, vp)
            implicit none
            real*8, intent(out)              :: a, e, s, o, vp
            real*8, dimension(5), intent(in) :: yi
            
            a  = yi(1) 
            e  = yi(2) 
            s  = yi(3) 
            o  = yi(4) 
            vp = yi(5)
        end subroutine get_from_yi
        
        subroutine get_from_big_y (y, y0, y1, y2)
            implicit none
            real*8, dimension(12), intent(in) :: y
            real*8, dimension(2), intent(out) :: y0
            real*8, dimension(5), intent(out) :: y1, y2
            
            y0 = y(1: 2)
            y1 = y(3: 7)
            y2 = y(8:12)
        end subroutine get_from_big_y
            
        subroutine get_from_ytidal (ytidal, y0, y1)
            implicit none
            real*8, dimension(7), intent(in)  :: ytidal
            real*8, dimension(2), intent(out) :: y0
            real*8, dimension(5), intent(out) :: y1
                       
            y0 = ytidal(1:2)
            y1 = ytidal(3:7)
        end subroutine get_from_ytidal
    
        
    ! DERIVATE

        subroutine dydtidal (y0, y1, planet, y01t, y10t)
            implicit none
            real*8, dimension(2), intent(in)          :: y0
            real*8, dimension(5), intent(in)          :: y1
            integer, intent(in)                       :: planet
            real*8                                    :: a1, e1, s1, o1, vp1, s0, o0
            real*8                                    :: fe1, fe2, fe3, fe4, fe5
            real*8                                    :: cos0, sin0, cos1, sin1
            real*8                                    :: n1
            real*8, dimension(size (y0)), intent(out) :: y01t ! d y0 / dt
            real*8, dimension(size (y1)), intent(out) :: y10t ! d y1 / dt
            integer                                   :: i
            
            i = planet + 1

            call get_from_y0 (y0, s0, o0)
            call get_from_yi (y1, a1, e1, s1, o1, vp1)
            
            n1 = n_f (a1, mu(i))
            call set_fl (e1, fe1, fe2, fe3, fe4, fe5) ! Can be "set_fs" [short] if desired
            call set_cosino (o0, cos0, sin0)
            call set_cosino (o1, cos1, sin1)

            y01t(1) = - TSk(1,i) * (fe1 * 0.5 * (1. + cos0**2) * s0 - fe2 * cos0 * n1) / &
                     & (C(1) * a1**6)            ! s0
            y01t(2) = TSk(1,i) * (fe1 * cos0 * 0.5 - fe2 * n1 / s0) * &
                     & sin0 / C(1) / a1**6       ! o0
            
            y10t(1) = (TSk(1,i) * (fe2 * cos0 * s0 / n1 - fe3) + &
                     & TSk(i,1) * (fe2 * cos1 * s1 / n1 - fe3)) * &
                      & 2. / mp(i) / a1**7       !a1
            y10t(2) = (TSk(1,i) * (11./18. * fe4 * cos0 * s0 / n1 - fe5) + &
                     & TSk(i,1) * (11./18. * fe4 * cos1 * s1 / n1 - fe5)) * &
                      & 9. * e1 / mp(i) / a1**8  !e1
            y10t(3) = -TSk(i,1) * (fe1 * (1. + cos1**2) * s1 * 0.5 - fe2 * cos1 * n1) / &
                       & C(i) / a1**6            !s1
            y10t(4) = TSk(i,1) * (fe1 * cos1 * 0.5 - fe2 * n1 / s1) * &
                      & sin1 / C(i) / a1**6      !o1   
            y10t(5) = 0. ! vp1
        end subroutine dydtidal

        function dydtidall (t, y) result (ynew)
            implicit none
            real*8, intent(in)               :: t
            real*8, dimension(:), intent(in) :: y
            real*8, dimension(size (y))      :: ynew
            real*8, dimension(2)             :: y0
            real*8, dimension(5)             :: y1, y2
            real*8, dimension(size (y0))     :: y01t, y02t
            real*8, dimension(size (y1))     :: y10t, y20t, y12g, y21g

            ynew = 0.

            call get_from_big_y (y, y0, y1, y2)
            call dydtidal (y0, y1, 1, y01t, y10t)
            ! call dydtidal (y0, y2, 2, y02t, y20t)
            ! call dydtgrav (mu1, mu2, m1, m2, y1, y2, y12g, y21g)

            ynew(1: 2) = y01t !+ y02t 
            ynew(3: 7) = y10t !+ y12g
            ! ynew(8:12) = y20t !+ y21g
        end function dydtidall

end module tidall_ip

! ! dydtidal0 (y0__, y1__) = ynew__
! subroutine dydtidal (m0, k0, c0, m1, k1, c1, mu1, m1p, y0, y1, y01t, y10t)
!     implicit none
!     real*8, intent(in)                        :: mu1, m0, k0, c0, m1, k1, c1, m1p
!     real*8, dimension(2), intent(in)          :: y0
!     real*8, dimension(5), intent(in)          :: y1
!     real*8                                    :: a1, e1, s1, o1, w1, s0, o0
!     real*8                                    :: fe1, fe2, fe3, fe4, fe5
!     real*8                                    :: cos0, sin0, cos1, sin1
!     real*8                                    :: n1, AM
!     real*8, dimension(size (y0)), intent(out) :: y01t ! d y0 / dt
!     real*8, dimension(size (y1)), intent(out) :: y10t ! d y1 / dt
    
!     call get_from_y0 (y0, s0, o0)
!     call get_from_yi (y1, a1, e1, s1, o1, w1)
    
!     n1 = n_f (a1, mu1)
!     AM = AngMom_f (a1, e1, mu1)
!     call set_fl (e1, fe1, fe2, fe3, fe4, fe5) ! Can be "set_fs" [short] if desired
!     call set_cosino (o0, cos0, sin0)
!     call set_cosino (o1, cos1, sin1)

!     y01t(1) = - k0 * n1 / (c0 * a1**6) * & ! s0
!             & (fe1 * 0.5 * (1. + cos0**2) * s0 / n1 - &
!             & fe2 * cos0)
!     y01t(2) = k0 * n1 / (c0 * s0 * a1**6) * & ! o0
!             & (fe1 * 0.5 * s0 / n1 * (cos0 - (c0 * s0) / (m1 * AM)) -  &
!             & fe2) * sin0

!     y10t(1) = 2. / m1p * a1**(-7) * & ! a1
!             & ((fe2 / n1) * (k0 * cos0 * s0 + k1 * cos1 * s1) - &
!             & fe3 * (k0 + k1))
!     y10t(2) = 9. / m1p * a1**(-8) * e1 * & ! e1
!             & ((11. / 18. * fe4 / n1) * (k0 * cos0 * s0 + k1 * cos1 * s1) - &
!             & fe5 * (k0 + k1))
!     y10t(3) = - k1 * n1 / (c1 * a1**6) * & !s1
!             & (fe1 * 0.5 * (1. + cos1**2) * s1 / n1 - &
!             & fe2 * cos1)
!     y10t(4) = k1 * n1 / (c1 * s1 * a1**6) * & !o1
!             & (fe1 * 0.5 * s1 / n1 * (cos1 - (c1 * s1) / (m0 * AM)) -  &
!             & fe2) * sin1
!     y10t(5) = 0. ! vp1
! end subroutine dydtidal


! ! dydtidal0 (y0__, y1__) = ynew__
! subroutine dydtgrav (mu1, mu2, m1, m2, y1, y2, y12g, y21g)
!     implicit none
!     real*8, intent(in)                        :: m1, m2, mu1, mu2
!     real*8, dimension(5), intent(in)          :: y1, y2
!     real*8                                    :: a1, e1, s1, o1, w1
!     real*8                                    :: a2, e2, s2, o2, w2
!     real*8                                    :: cos1, sin1, cos2, sin2
!     real*8                                    :: cte1, cte2
!     real*8                                    :: aux
!     real*8, dimension(size (y1)), intent(out) :: y12g, y21g ! d y1 / dt ; d y2 / dt
    
!     call get_from_yi (y1, a1, e1, s1, o1, w1)
!     call get_from_yi (y2, a2, e2, s2, o2, w2)

!     call set_cosino (o1, cos1, sin1)
!     call set_cosino (o2, cos2, sin2)
    
!     cte1 = m2 / (a1**2 * e1 * n_f (a1, mu1))
!     cte2 = m1 / (a2**2 * e2 * n_f (a2, mu2))

!     aux = 3. * a1**2 * G / (16. * a2**4 * (1. - e2**2)**(2.5))

!     y12g(1) = 0.
!     y12g(2) = aux * 5. * a1 * e1 * e2 * sin (w1 - w2) * cte1
!     y12g(3) = 0.
!     y12g(4) = 0.
!     y12g(5) = aux * (4. * a2 * e1 * (1. - e2**2) - 5 * a1 * e2 * cos (w1 - w2)) * cte1

!     y21g(1) = 0.
!     y21g(2) = - y12g(2) * cte2 / cte1
!     y21g(3) = 0.
!     y21g(4) = 0.
!     y21g(5) = aux * (2. * a2 * e2 * (2. + 3. * e1**2 + e2**2) - &
!                    & (5. * a1 * e1 * (1. + 4. * e2**2) * cos (w1 - w2)) / (1. - e2**2)) * cte2
! end subroutine dydtgrav

! function dydtidall (t, y) result (ynew)
!     implicit none
!     real*8, intent(in)               :: t
!     real*8, dimension(:), intent(in) :: y
!     real*8, dimension(size (y))      :: ynew
!     real*8, dimension(2)             :: y0
!     real*8, dimension(5)             :: y1, y2
!     real*8, dimension(size (y0))     :: y01t, y02t
!     real*8, dimension(size (y1))     :: y10t, y20t, y12g, y21g

!     call get_from_big_y (y, y0, y1, y2)
!     call dydtidal (m(1), TSk(1, 2), C(1), &
!                  & m(2), TSk(2, 1), C(2), mu(2), mp(2), &
!                  & y0, y1, y01t, y10t)
!     call dydtidal (m(1), TSk(1, 3), C(1), &
!                  & m(3), TSk(3, 1), C(3), mu(3), mp(3), &
!                  & y0, y2, y02t, y20t)
!     !call dydtgrav (mu1, mu2, m1, m2, y1, y2, y12g, y21g)

!     ynew(1: 2) = y01t + y02t 
!     ynew(3: 7) = y10t !+ y12g
!     ynew(8:12) = y20t !+ y21g
! end function dydtidall
