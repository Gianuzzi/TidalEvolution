program tidal

use run
use tidall
use integrators
use bstoer

implicit none
!------------------ INPUT PARAMETERS ------------------
! Objects
!! Object 0
m(1)   = 1.0             ! [Ms]
s0     = TWOPI / 28.     ! [rad day⁻¹]
o0     = 25. * DEG2RAD   ! [rad]
r(1)   = RS2UA           ! [UA]
alp(1) = 0.4             ! []
Q(1)   = 1.e6            ! [?]

!! Object 1
m(2)   = MT2MS           ! [Ms]
a1     = 0.05            ! [UA] 
e1     = 0.1             ! []
s1     = TWOPI / 0.01    ! [rad day⁻¹]
o1     = 80. * DEG2RAD   ! [rad]
vp1    = 0.              ! [rad]
r(2)   = RT2UA           ! [UA]
alp(2) = 0.4             ! []
Q(2)   = 1.e2            ! [?]

!! Object 2
m(3)   = MJ2MS           ! [Ms]
a2     = 0.2             ! [UA] 
e2     = 0.1             ! []
s2     = TWOPI / 0.01    ! [rad day⁻¹]
o2     = 40. * DEG2RAD   ! [rad]
vp2    = PI              ! [rad]
r(3)   = RJ2UA           ! [UA]
alp(3) = 0.4             ! []
Q(3)   = 1.e5            ! [?]

! Run conditions
t0       = 0. * YR2DAY    ! [days]
dt       = 0.0001 * YR2DAY  ! [days] ![First & min]
tf       = 5.e10 * YR2DAY ! [days]
n_points = 5000           ! N_output

! Integration conditions
beta   = 0.95   ! Learning rate
e_tol  = 1e-12  ! Approx Absolute e_calc (|Ysol - Ypred|)

! Output
filename = "Salida4.txt"
!------------------------------------------------------

!-------------- SET DERIVED PARAMETERS --------------
! Set almos everything
call set_initial_params ((/a1, a2/), Q, alp, r, m, C, mu, mp, n, k2dt, TSk)

! Calculated LOOP parameters
n_iter = int ((tf - t0) / dt, kind=8)
Logt   = exp (log (tf - t0) / (n_points - 1))
t_add  = Logt
print *, "Approximate Iterations:", n_iter
!------------------------------------------------------

!------------------------ FILE ------------------------
inquire (file=trim (filename), exist=file_exists)

if (file_exists) then
    print '("File ",A, " already exists.")', trim (filename)
    print*, "Overwrite file:"
    print*, "    [Y]: Yes."
    print*, "      N: No."
    read*, selection
    select case (selection)
        case ("Y")
            print*, "Overwriting..."
        case ("y")
            print*, "Overwriting..."
        case default
            print*, ("Exiting.")
            stop (0)
    end select
end if

open (20, file=trim (filename)//".params", status='replace')
write (20, '(17(A,1X))') "m0", "r0", "al0", "Q0", &
                       & "m1", "r1", "al1", "Q1", &
                       & "m2", "r2", "al2", "Q2", &
                       & "s0", "o0", &
                       & "a1", "e1", "s1", "o1", "vp1", &
                       & "a2", "e2", "s2", "o2", "vp2", &
                       & "t0", "dt", "tf"
write (20,'(17(E16.9,1X))') m(1), r(1), alp(1), Q(1), &
                          & m(2), r(2), alp(2), Q(2), &
                          & m(3), r(3), alp(3), Q(3), &
                          & s0, o0, &
                          & a1, e1, s1, o1, vp1, &
                          & a2, e2, s2, o2, vp2, &
                          & t0, dt, tf
close (20)
!------------------------------------------------------

!------------------------ START -----------------------

! Time
call cpu_time (start_time)

! Output file
open (10, file=trim (filename), status='replace')
!------------------------------------------------------

!------------------------ LOOP ------------------------

! Initial parameters
call set_y0 (s0, o0, y0)
call set_yi (a1, e1, s1, o1, vp1, y1)
call set_yi (a2, e2, s2, o2, vp2, y2)
call set_big_y (y0, y1, y2, y)
t       = t0
t_out   = t0 + t_add
dt_min  = dt
dt_adap = dt ! For adaptive step
i       = 0

! Do loop
do while (t < tf)
    
    ! CHECK
    !! s0, o0, a1, K1, s1, o1, H1, a2, K2, s2, o2, H2
    !!  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12
    if (min (y(3), y(8)) < r(1) * 4.) then
        write (*,*) "End of RUN. [Encounter]"
        exit
    end if

    ! Output
    if (t >= t_out) then
        t_add = t_add * Logt
        t_out = t0 + t_add
        write (*, "(I11, 4(E14.4, 1X))") i, t / YR2DAY, dt / YR2DAY
        ! write (*, "(A8, 12(E14.4, 1X))") "Valores:", y
        write (10,*) y, n_f (y(3), mu(2)), n_f (y(8), mu(3)), t, dt
    end if
    
    !!! Execute an integration method (uncomment one of theese)
    !  call integ_caller (t, y, dt, dydtidall, rungek4, ynew)
    !  call rec_rk_adap (t, y, dt_adap, dydtidall, rungek4, 4, e_tol, beta, dt_min, dt, ynew)

      call Verner5_6 (t, y, dt_adap, dydtidall, e_tol, beta, dt_min, dt, ynew)
    
    !! Modulate and avoid too small angles
    ynew(2)  = max (1.0d-15, mod (ynew(2), TWOPI))
    ynew(6)  = max (1.0d-15, mod (ynew(6), TWOPI))
    ynew(11) = max (1.0d-15, mod (ynew(11), TWOPI))
    
    ! Update parameters
    i = i + 1
    t = t + dt
    y = ynew
    
end do
!------------------------------------------------------

!------------------------ END -------------------------

! Output file
close (10)

! Time
call cpu_time (final_time)

! Messages
print '("Total iterations:",I11)', i
print '("Final time:",E16.6," [yrs]")', t / YR2DAY
print '("Total running time:",F10.4," [sec].")', final_time - start_time
!------------------------------------------------------

end program tidal


