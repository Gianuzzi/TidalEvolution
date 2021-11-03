module bstoer

use tidall, only: dydtidall

contains

subroutine derivs (t, y, f)
implicit real*8 (a-h,k-z)
real*8 y(12),f(12)

f = dydtidall (t, y)

end subroutine derivs

subroutine bs (t0,y,delt,in,eps,step,ynew)
implicit real*8 (a-h,k-z)
real*8 y(12),f(12),ynew(12), time
!
ynew = y
time = t0
tfinal = time + delt
htry   = step
do while (time < tfinal)
    if (time+htry > tfinal) htry = tfinal - time
    call derivs (time,ynew,f)
    call bstep  (ynew,f,in,time,htry,eps,hdid,hnext)
    htry = hnext
end do
step = hnext
!
end subroutine bs


subroutine bstep (y,dydx,inv,x,htry,eps,hdid,hnext)
implicit real*8 (a-h,o-z)
parameter (imax=12,nmax=12,kmaxx=8,safe1=.25,safe2=.7, &
    redmax=1.e-5,redmin=.7,tini=1.e-30,scalmx=.1)
integer nseq(kmaxx+1)
real*8 y(nmax),dydx(nmax),yscal(nmax),a(kmaxx+1),alf(kmaxx,kmaxx)
real*8 err(kmaxx),yerr(nmax),ysav(nmax),yseq(nmax)
logical first,reduct
save a,alf,epsold,first,kmax,kopt,nseq,xnew
data first/.true./,epsold/-1./
data nseq /2,4,6,8,10,12,14,16,18/
!
if (eps.ne.epsold) then
    hnext = -1.0d29
    xnew  = -1.0d29
    eps1  = safe1*eps
    a(1) = nseq(1) + 1
    do k = 1,kmaxx
    a(k+1) = a(k) + nseq(k+1)
    end do
    do iq = 2,kmaxx
    do k = 1,iq-1
        alf(k,iq)=eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.)*(2*k+1)))
    end do
    end do
    epsold = eps
    do kopt = 2,kmaxx-1
    if (a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt)) goto 1
    end do
1   kmax = kopt
end if
h = htry
do i = 1,inv
    ysav(i) = y(i)
end do
if (h.ne.hnext.or.x.ne.xnew) then
    first = .true.
    kopt = kmax
end if
reduct = .false.
2 do k = 1,kmax
    xnew = x + h
    if (xnew.eq.x) return
    call mmid (ysav,dydx,inv,x,h,nseq(k),yseq)
    xest = (h/nseq(k))**2
    call pzextr (k,xest,yseq,y,yerr,inv)
    if (k.ne.1) then
    errmax = tini
    do i = 1,inv
        errmax = max(errmax,abs(yerr(i)))
    end do
    errmax = errmax/eps
    km = k - 1
    err(km) = (errmax/safe1)**(1./(2*km+1))
    end if
    if (k.ne.1.and.(k.ge.kopt-1.or.first)) then
    if (errmax.lt.1.) goto 4
    if (k.eq.kmax.or.k.eq.kopt+1) then
        red = safe2/err(km)
        goto 3
    else if (k.eq.kopt) then
        if (alf(kopt-1,kopt).lt.err(km)) then
        red = 1./err(km)
        goto 3
        end if
    else if (kopt.eq.kmax) then
        if (alf(km,kmax-1).lt.err(km)) then
        red = alf(km,kmax-1)*safe2/err(km)
        goto 3
        end if
    else if (alf(km,kopt).lt.err(km)) then
        red = alf(km,kopt-1)/err(km)
        goto 3
    end if
    end if
end do
3 red = min(red,redmin)
red = max(red,redmax)
h = h*red
reduct = .true.
goto 2
4 x = xnew
hdid = h
first = .false.
wrkmin = 1.d35
do kk = 1,km
    fact = max(err(kk),scalmx)
    work = fact*a(kk+1)
    if (work.lt.wrkmin) then
    scala = fact
    wrkmin = work
    kopt = kk + 1
    end if
end do
hnext = h/scala
if (kopt.ge.k.and.kopt.ne.kmax.and..not.reduct) then
    fact = max(scala/alf(kopt-1,kopt),scalmx)
    if (a(kopt+1)*fact.le.wrkmin) then
    hnext = h/fact
    kopt = kopt + 1
    end if
end if
!
end subroutine bstep


subroutine mmid (y,dydx,nvar,xs,htot,nstep,yout)
implicit real*8 (a-h,o-z)
parameter (imax=12,nmax=12)
real*8 dydx(nmax),y(nmax),ym(nmax),yn(nmax),yout(nmax),ders(nmax)
!
h  = htot/nstep
do i = 1,nvar
    ym(i) = y(i)
    yn(i) = y(i) + h*dydx(i)
end do
x  = xs + h
call derivs (x,yn,ders)
h2 = 2.*h
do n = 2,nstep
    do i = 1,nvar
    swap = ym(i) + h2*ders(i)
    ym(i) = yn(i)
    yn(i) = swap
    end do
    x = x + h
    call derivs (x,yn,ders)
end do
do i = 1,nvar
    yout(i) = 0.5*(ym(i)+yn(i)+h*ders(i))
end do
!    
end subroutine mmid


subroutine pzextr (iest,xest,yest,yz,dy,nv)
implicit real*8 (a-h,o-z)
parameter (imax=12,nmax=12,immax=23)
real*8 dy(nmax),yest(nmax),yz(nmax),d(nmax),qcol(nmax,immax)
real*8 x(immax)
save qcol,x
!    
x(iest) = xest
do j = 1,nv
    dy(j) = yest(j)
    yz(j) = yest(j)
end do
if (iest.eq.1) then
    do j = 1,nv
    qcol(j,1) = yest(j)
    end do
else
    do j = 1,nv
    d(j) = yest(j)
    end do
    do k1 = 1,iest-1
    delta = 1.0/(x(iest-k1)-xest)
    f1 = xest*delta
    f2 = x(iest-k1)*delta
    do j = 1,nv
        q = qcol(j,k1)
        qcol(j,k1) = dy(j)
        delta = d(j) - q
        dy(j) = f1*delta
        d(j)  = f2*delta
        yz(j) = yz(j) + dy(j)
    end do
    end do
    do j = 1,nv
    qcol(j,iest) = dy(j)
    end do
end if
! 
end subroutine pzextr

end module bstoer
