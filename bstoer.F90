module bstoer

    implicit none

    abstract interface
        ! f__ (t, y__) = (f_i (t, y__), ..., f_n (t, y__)) = der__
        function dydt_tem (t, y) result (der)
            implicit none
            real*8, intent(in)               :: t
            real*8, dimension(:), intent(in) :: y
            real*8, dimension(size (y))      :: der
            ! Here must be every f_i defined explicitly
        end function dydt_tem
    end interface

    contains
    
        !---------------------------------
        !   Call Bulirsch_Stoer
        !---------------------------------
        
        subroutine BStoer_caller (t, y, dt_adap, dydt, e_tol, dt_min, dt, ynew)
            implicit none
            real*8, intent(in)                       :: t, e_tol, dt_min, dt
            real*8, dimension(:), intent(in)         :: y
            real*8, intent(inout)                    :: dt_adap
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: yaux
            real*8                                   :: time, t_end, dtmin, dtused

            ynew  = y
            time  = t
            t_end = time + dt
            dtmin = min (dt_min, dt)
            do while (time < t_end)
                yaux  = ynew
                dt_adap = min (dt_adap, t_end - time)
                call Bulirsch_Stoer (time, yaux, dt_adap, dydt, e_tol, dtused, ynew)
                time = time + dtused
            end do
        end subroutine BStoer_caller

        !! Bulirsch Stoer main integrator
    
        subroutine Bulirsch_Stoer (t, y, dt_adap, dydt, e_tol, dt_used, ynew)
            implicit none
            real*8, intent(in)                       :: t, e_tol
            real*8, intent(inout)                    :: dt_adap, dt_used
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8                                   :: time, dtry
            integer*4, save                          :: sizey
            
            sizey = size (y)
            time  = t
            dtry  = dt_adap
            call bstep (ynew, dydt, sizey, time, dtry, e_tol, dt_used, dt_adap)
        end subroutine Bulirsch_Stoer

        !!

        !! Auxiliar subroutines

        subroutine bstep (y, dydt, sizey, x, htry, eps, hdid, hnext)
            implicit none        
            procedure(dydt_tem)   :: dydt
            integer*4, intent(in) :: sizey
            real*8, intent(in)    :: htry, eps
            real*8, intent(inout) :: x, hdid, hnext
            real*8, dimension(sizey), intent(inout) :: y
            real*8, dimension(sizey)                :: yerr, ysav, yseq, der
            real*8, dimension(kmaxx * 2)        :: xpz
            real*8, dimension(sizey, kmaxx * 2) :: qcolpz
            integer*4, parameter :: kmaxx = 8
            real*8, parameter    :: safe1 = .25, safe2 = .7
            real*8, parameter    :: redmax = 1.e-5, redmin = .7
            real*8, parameter    :: tini = 1.e-30, scalmx = .1
            real*8, dimension(kmaxx)                   :: err
            real*8, dimension(kmaxx + 1), save         :: arr
            real*8, dimension(kmaxx, kmaxx), save      :: alf
            integer*4, parameter, dimension(kmaxx + 1) :: nseq = (/2, 4, 6, 8, 10, 12, 14, 16, 18/)
            logical         :: reduct
            logical, save   :: first = .true.
            integer*4, save :: kmax, kopt
            real*8, save    :: xnew, epsold = -1.
            real*8          :: wrkmin, fact, work, scala, eps1, xest, errmax, red, h
            integer         :: k, iq, i ,km, kk

            der = dydt (x, y)

            if (eps .ne. epsold) then
                hnext = -1.0e29
                xnew  = -1.0e29
                eps1  = safe1 * eps
                arr(1) = nseq(1) + 1
                do k = 1, kmaxx
                    arr(k + 1) = arr(k) + nseq(k + 1)
                end do
                do iq = 2, kmaxx
                    do k = 1, iq - 1
                        alf(k, iq) = eps1**((arr(k + 1) - arr(iq + 1)) / ((arr(iq + 1) - arr(1) + 1.) * (2. * k + 1)))
                    end do
                end do
                epsold = eps
                do kopt = 2, kmaxx - 1
                    if (arr(kopt + 1) .gt. arr(kopt) * alf(kopt - 1, kopt)) then
                        goto 1
                    end if
                end do
                1 kmax = kopt
            end if
            h = htry
            do i = 1, sizey
                ysav(i) = y(i)
            end do
            if (h .ne. hnext .or. x .ne. xnew) then
                first = .true.
                kopt = kmax
            end if
            reduct = .false.
            2 do k = 1, kmax
                xnew = x + h
                if (xnew .eq. x) then
                    return
                end if
                call mmid (ysav, dydt, der, sizey, x, h, nseq(k), yseq)
                xest = (h / nseq(k))**2
                call pzextr (k, xest, yseq, y, yerr, sizey, qcolpz, xpz)
                if (k .ne. 1) then
                    errmax = tini
                    do i = 1, sizey
                        errmax = max (errmax, abs (yerr(i)))
                    end do
                    errmax = errmax / eps
                    km = k - 1
                    err(km) = (errmax / safe1)**(1. / (2. * km + 1))
                end if
                if (k .ne. 1 .and. (k .ge. kopt - 1 .or. first)) then
                    if (errmax .lt. 1.) then
                        goto 4
                    end if
                    if (k .eq. kmax .or. k .eq. kopt + 1) then
                        red = safe2 / err(km)
                        goto 3
                    else if (k .eq. kopt) then
                        if (alf(kopt - 1, kopt) .lt. err(km)) then
                            red = 1. / err(km)
                            goto 3
                        end if
                    else if (kopt .eq. kmax) then
                        if (alf(km, kmax - 1) .lt. err(km)) then
                            red = alf(km, kmax - 1) * safe2 / err(km)
                            goto 3
                        end if
                    else if (alf(km, kopt) .lt. err(km)) then
                        red = alf(km, kopt - 1) / err(km)
                        goto 3
                    end if
                end if
            end do
            3 red = min (red, redmin)
            red = max (red, redmax)
            h = h * red
            reduct = .true.
            goto 2
            4 x = xnew
            hdid = h
            first = .false.
            wrkmin = 1.e35
            do kk = 1, km
                fact = max (err(kk), scalmx)
                work = fact * arr(kk + 1)
                if (work .lt. wrkmin) then
                    scala = fact
                    wrkmin = work
                    kopt = kk + 1
                end if
            end do
            hnext = h / scala
            if (kopt .ge. k .and. kopt .ne. kmax .and. .not. reduct) then
                fact = max (scala / alf(kopt - 1, kopt), scalmx)
                if (arr(kopt + 1) * fact .le. wrkmin) then
                    hnext = h / fact
                    kopt = kopt + 1
                end if
            end if
        end subroutine bstep

        subroutine mmid (y, dydt, dydx, sizey, xs, htot, nstep, yout)
            integer*4, intent(in) :: sizey, nstep
            procedure(dydt_tem)   :: dydt
            real*8, intent(in)    :: xs, htot
            real*8, dimension(sizey)              :: ym, yn, dernew
            real*8, dimension(sizey), intent(in)  :: y, dydx
            real*8, dimension(sizey), intent(out) :: yout
            real*8  :: x, h, h2, swap
            integer :: i, n 

            h  = htot / real (nstep, kind=8)
            do i = 1, sizey
                ym(i) = y(i)
                yn(i) = y(i) + h * dydx(i)
            end do
            x  = xs + h
            dernew = dydt(x, yn)
            h2 = 2. * h
            do n = 2, nstep
                do i = 1, sizey
                    swap = ym(i) + h2 * dernew(i)
                    ym(i) = yn(i)
                    yn(i) = swap
                end do
                x = x + h
                dernew = dydt(x, yn)
            end do
            do i = 1, sizey
                yout(i) = 0.5 * (ym(i) + yn(i) + h * dernew(i))
            end do
        end subroutine mmid

        subroutine pzextr (iest, xest, yest, yz, dy, sizey, qcol, x)
            integer*4, intent(in) :: iest, sizey
            real*8, intent(in)    :: xest
            real*8, dimension(sizey)                 :: d
            real*8, dimension(sizey), intent(in)     :: yest
            real*8, dimension(sizey), intent(inout)  :: dy, yz
            integer*4, parameter                     :: IEST_MAX=16
            real*8, dimension(IEST_MAX), intent(inout)        :: x
            real*8, dimension(sizey, IEST_MAX), intent(inout) :: qcol
            integer :: j, k1
            real*8  :: delta, f1, f2, q

            x(iest) = xest
            do j = 1, sizey
                dy(j) = yest(j)
                yz(j) = yest(j)
            end do
            if (iest .eq. 1) then
                do j = 1, sizey
                    qcol(j, 1) = yest(j)
                end do
            else
                do j = 1, sizey
                    d(j) = yest(j)
                end do
                do k1 = 1, iest - 1
                    delta = 1.0 / (x(iest - k1) - xest)
                    f1 = xest * delta
                    f2 = x(iest - k1) * delta
                    do j = 1, sizey
                        q = qcol(j, k1)
                        qcol(j, k1) = dy(j)
                        delta = d(j) - q
                        dy(j) = f1 * delta
                        d(j)  = f2 * delta
                        yz(j) = yz(j) + dy(j)
                    end do  
                end do
                do j = 1, sizey
                    qcol(j, iest) = dy(j)
                end do
            end if
        end subroutine pzextr

        !!

end module bstoer
