module tidall

    use const

    implicit none
    integer, parameter :: N_bodies = 3, y0_s = 2, yi_s = 5
    integer, parameter :: y_s = y0_s + yi_s * (N_bodies - 1)

    real*8, dimension(N_bodies) :: m            ! Masses (constants)
    real*8, dimension(N_bodies) :: r, alp, C    ! Bodies constant parameters
    real*8, dimension(N_bodies) :: Q, k2dt      ! Bodies derived parameters
    real*8, dimension(N_bodies) :: mu, mp, n    ! Bodies 1... derived parameters
    
    real*8, dimension(N_bodies, N_bodies) :: TSk ! TIDAL STRENGTH [k]

    real*8 :: s0, o0                         ! Body 0 parameters
    real*8 :: a1, e1, s1, o1, vp1            ! Body 1 parameters
    real*8 :: a2, e2, s2, o2, vp2            ! Body 2 parameters

    real*8, dimension(y0_s) :: y0       ! Body 0 Parameters array
    real*8, dimension(yi_s) :: y1, y2   ! Body 1 and Body 2 Parameters array
    real*8, dimension(y_s)  :: y, ynew  ! Parameters array, and same for iteration
    
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
    
    ! SETTERS
        subroutine set_initial_params (a, Q, alpha, radius, m, C, mu, mp, n, k2dt, K)
            implicit none
            real*8, dimension(N_bodies - 1), intent(in)        :: a
            real*8, dimension(N_bodies), intent(in)            :: m, Q, alpha, radius
            real*8, dimension(N_bodies), intent(out)           :: mu, mp, n, C, k2dt
            real*8, dimension(N_bodies, N_bodies), intent(out) :: K
            integer                                            :: i, j
            
            K = 0.
            do i = N_bodies, 1, -1 !Backwards, for n(1) setting
                if (i > 1) then
                    mp(i) = mprime_f (sum (m(:(i - 1))), m(i))
                    mu(i) = mu_f (sum (m(:(i - 1))), m(i))                    
                    n(i)  = n_f (a(i - 1), mu(i))
                else
                    mu(i) = 0.   ! Unused
                    mp(i) = 0.   ! Unused
                    n(i)  = n(2) ! n(1) = n(2)
                    
                end if
                C(i)    = alpha(i) * m(i) * radius(i)**2
                k2dt(i) = 1.5 / Q(i) / n(i)                
                do j = 1, N_bodies
                    if (i .ne. j) then
                        K(i, j) = TidalStrength_f (n(i), m(j), radius(i), Q(i))
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
            real*8, intent(in)                   :: s0, o0
            real*8, dimension(y0_s), intent(out) :: y0
            
            y0(1) = s0
            y0(2) = o0
        end subroutine set_y0

        subroutine set_yi (a, e, s, o, vp, yi)
            implicit none
            real*8, intent(in)                   :: a, e, s, o, vp
            real*8                               :: K, H
            real*8, dimension(yi_s), intent(out) :: yi
            
            call get_KH (e, vp, K, H)
            yi(1) = a
            yi(2) = K
            yi(3) = s
            yi(4) = o
            yi(5) = H
        end subroutine set_yi
        
        subroutine set_big_y (y0, y1, y2, y)
            implicit none
            real*8, dimension(y0_s), intent(in) :: y0
            real*8, dimension(yi_s), intent(in) :: y1, y2
            real*8, dimension(y_s), intent(out) :: y
            
            y(: y0_s)                 = y0
            y(y0_s + 1 : y0_s + yi_s) = y1
            y(y0_s + yi_s + 1 :)      = y2
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

        subroutine get_KH(e, vp, K, H)
            implicit none
            real*8, intent(in)  :: e, vp
            real*8, intent(out) :: K, H
            
            K = e * cos (vp)
            H = e * sin (vp)
        end subroutine get_KH

        subroutine get_evp(K, H, e, vp)
            implicit none
            real*8, intent(in)  :: K, H
            real*8, intent(out) :: e, vp

            e  = sqrt (K**2 + H**2)
            vp = atan2 (H, K)        
        end subroutine get_evp
            
        subroutine get_from_y0 (y0, s0, o0)
            implicit none
            real*8, intent(out)                 :: s0, o0
            real*8, dimension(y0_s), intent(in) :: y0
            
            s0 = y0(1)
            o0 = y0(2)
        end subroutine get_from_y0
    
        subroutine get_from_yi (yi, a, K, s, o, H)
            implicit none
            real*8, dimension(yi_s), intent(in) :: yi
            real*8, intent(out)                 :: a, K, s, o, H
            
            a = yi(1) 
            K = yi(2) 
            s = yi(3) 
            o = yi(4) 
            H = yi(5)
        end subroutine get_from_yi
        
        subroutine get_from_big_y (y, y0, y1, y2)
            implicit none
            real*8, dimension(y_s), intent(in)   :: y
            real*8, dimension(y0_s), intent(out) :: y0
            real*8, dimension(yi_s), intent(out) :: y1, y2
            
            y0 = y(: y0_s)
            y1 = y(y0_s + 1 : y0_s + yi_s)
            y2 = y(y0_s + yi_s + 1 :)
        end subroutine get_from_big_y
            
        subroutine get_from_ytidal (ytidal, y0, y1)
            implicit none
            real*8, dimension(y0_s + yi_s), intent(in) :: ytidal
            real*8, dimension(y0_s), intent(out)       :: y0
            real*8, dimension(yi_s), intent(out)       :: y1
                       
            y0 = ytidal(: y0_s)
            y1 = ytidal(y0_s + 1 :)
        end subroutine get_from_ytidal
    
        
    ! DERIVATE

        subroutine dydtidal (y0, y1, planet, y01t, y10t)
            implicit none
            real*8, dimension(y0_s), intent(in)  :: y0
            real*8, dimension(yi_s), intent(in)  :: y1
            integer, intent(in)                  :: planet
            real*8                               :: s0, o0
            real*8                               :: a1, K1, s1, o1, H1
            real*8                               :: e1, vp1
            real*8                               :: fe1, fe2, fe3, fe4, fe5
            real*8                               :: cos0, sin0, cos1, sin1
            real*8                               :: n1, de1dt
            real*8, dimension(y0_s), intent(out) :: y01t ! d y0 / dt
            real*8, dimension(yi_s), intent(out) :: y10t ! d y1 / dt
            integer                              :: i
            
            i = planet + 1

            call get_from_y0 (y0, s0, o0)
            call get_from_yi (y1, a1, K1, s1, o1, H1)

            call get_evp (K1, H1, e1, vp1)

            call set_fl (e1, fe1, fe2, fe3, fe4, fe5) ! Can be "set_fs" [short] if desired
            call set_cosino (o0, cos0, sin0)
            call set_cosino (o1, cos1, sin1)

            n1 = n_f (a1, mu(i))

            de1dt = ((TSk(1,i) * (11./18. * fe4 * cos0 * s0 / n1 - fe5) + &
                    & TSk(i,1) * (11./18. * fe4 * cos1 * s1 / n1 - fe5)) * &
                    & 9. * e1 / mp(i) / a1**8)

            y01t(1) = - TSk(1,i) * (fe1 * 0.5 * (1. + cos0**2) * s0 - fe2 * cos0 * n1) / &
                     & (C(1) * a1**6)            ! d(s0) / dt
            y01t(2) = TSk(1,i) * (fe1 * cos0 * 0.5 - fe2 * n1 / s0) * &
                     & sin0 / C(1) / a1**6       ! d(o0) / dt
            
            y10t(1) = (TSk(1,i) * (fe2 * cos0 * s0 / n1 - fe3) + &
                     & TSk(i,1) * (fe2 * cos1 * s1 / n1 - fe3)) * &
                     & 2. / mp(i) / a1**7        ! d(a1) / dt
            y10t(2) = de1dt * cos (vp1)          ! d(K1) / dt = d(e1) / dt * cos(vp1)
            y10t(3) = - TSk(i,1) * (fe1 * 0.5 * (1. + cos1**2) * s1 - fe2 * cos1 * n1) / &
                     & C(i) / a1**6              ! d(s1) / dt
            y10t(4) = TSk(i,1) * (fe1 * cos1 * 0.5 - fe2 * n1 / s1) * &
                     & sin1 / C(i) / a1**6       ! d(o1) / dt
            y10t(5) = de1dt * sin (vp1)          ! d(H1) / dt = d(e1) / dt * sin(vp1)
        end subroutine dydtidal
        
        ! dydtidal0 (y0__, y1__) = ynew__
        subroutine dydtgrav (y1, planet1, y2, planet2, y12g, y21g)
            implicit none
            real*8, dimension(yi_s), intent(in)  :: y1, y2
            integer, intent(in)                  :: planet1, planet2
            real*8                               :: a1, K1, s1, o1, H1
            real*8                               :: a2, K2, s2, o2, H2
            real*8                               :: extra1, extra2
            real*8                               :: alpha, factor
            integer                              :: i1, i2
            real*8, dimension(yi_s), intent(out) :: y12g, y21g ! d y1 / dt ; d y2 / dt
            
            i1 = planet1 + 1
            i2 = planet2 + 1

            call get_from_yi (y1, a1, K1, s1, o1, H1)
            call get_from_yi (y2, a2, K2, s2, o2, H2)
            
            alpha  = a1 / a2
            factor = G / a2**3
            extra1 = m(i2) / n_f (a1, mu(i1))
            extra2 = m(i1) * alpha**2 / n_f (a2, mu(i2))
            
            y12g = 0.
            y21g = 0.
            
            y12g(2) = - factor * (0.75 * H1 - 0.9375 * alpha * H2) * extra1 ! d(K1) / dt
            y12g(5) = factor * (0.75 * K1 - 0.9375 * alpha * K2) * extra1   ! d(H1) / dt
            
            y21g(2) = - factor * (0.75 * H2 - 0.9375 * alpha * H1) * extra2 ! d(K2) / dt
            y21g(5) = factor * (0.75 * K2 - 0.9375 * alpha * K1) * extra2   ! d(H2) / dt
        end subroutine dydtgrav

        subroutine dydtrela (yi, planet, yir)
            implicit none
            real*8, dimension(yi_s), intent(in)  :: yi
            integer, intent(in)                  :: planet
            real*8                               :: a, K, s, o, H
            real*8                               :: e, vp
            real*8                               :: dvpdt
            real*8, dimension(yi_s), intent(out) :: yir ! d y1 / dt
            integer                              :: i
            
            i = planet + 1

            call get_from_yi (yi, a, K, s, o, H)
            call get_evp (K, H, e, vp)

            dvpdt = 3. * mu(i) * n_f (a, mu(i)) / (a * C_SPEED**2)

            yir = 0.         

            yir(2) = - e * sin (vp) * dvpdt  ! d(K) / dt = - e * sin(vp) * d(vp) / dt 
            yir(5) = e * cos (vp) * dvpdt    ! d(H) / dt = e * cos(vp) * d(vp) / dt 
        end subroutine dydtrela

        function dydtidall (t, y) result (ynew)
            implicit none
            real*8, intent(in)               :: t
            real*8, dimension(:), intent(in) :: y
            real*8, dimension(size (y))      :: ynew
            real*8, dimension(y0_s)          :: y0
            real*8, dimension(yi_s)          :: y1, y2
            real*8, dimension(y0_s)          :: y01t, y02t
            real*8, dimension(yi_s)          :: y10t, y20t, y12g, y21g, y1r, y2r

            call get_from_big_y (y, y0, y1, y2)

            ynew = 0.
            y10t = 0.
            y20t = 0.
            y12g = 0.
            y21g = 0.
            y1r  = 0.
            y2r  = 0.

            call dydtidal (y0, y1, 1, y01t, y10t)
            call dydtidal (y0, y2, 2, y02t, y20t)
            call dydtgrav (y1, 1, y2, 2, y12g, y21g)
            call dydtrela (y1, 1, y1r)
            call dydtrela (y2, 2, y2r)

            ynew(: y0_s)                 = y01t + y02t       ! d(y0)/dt
            ynew(y0_s + 1 : y0_s + yi_s) = y10t + y12g + y1r ! d(y1)/dt
            ynew(y0_s + yi_s + 1 :)      = y20t + y21g + y2r ! d(y2)/dt
        end function dydtidall
        
end module tidall
