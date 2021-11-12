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
t0       = 0. * YR2DAY   ! [days]
dt_min   = 1e-4 * YR2DAY ! [days] ![For fixed dt integrators]
tf       = 5.e9 * YR2DAY ! [days]
n_points = 5000          ! Approx N_output

! Integration conditions
beta   = 0.95   ! Learning rate
e_tol  = 1e-12  ! Approx Absolute e_calc (|Ysol - Ypred|)

! Output
filename = "3_bodies"
!------------------------------------------------------

!-------------- SET DERIVED PARAMETERS --------------
! Set almost everything
call set_initial_params ((/a1, a2/), Q, alp, r, m, C, mu, mp, n, k2dt, TSk)

! Calculated LOOP parameters
Logt = exp (log (tf - t0) / (real (n_points, kind=8) - 1.))
!------------------------------------------------------

!------------------------ FILE ------------------------
inquire (file=trim (filename)//".txt", exist=file_exists)

if (file_exists) then
    print '("File ",A, " already exists.")', trim (filename)//".txt"
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

open (20, file=trim (filename)//"_params.txt", status='replace')
write (20, '(27(A,1X))') "m0", "r0", "al0", "Q0", &
                       & "m1", "r1", "al1", "Q1", &
                       & "m2", "r2", "al2", "Q2", &
                       & "s0", "o0", &
                       & "a1", "e1", "s1", "o1", "vp1", &
                       & "a2", "e2", "s2", "o2", "vp2", &
                       & "t0", "dt_min", "tf"
write (20,'(27(E16.9,1X))') m(1), r(1), alp(1), Q(1), &
                          & m(2), r(2), alp(2), Q(2), &
                          & m(3), r(3), alp(3), Q(3), &
                          & s0, o0, &
                          & a1, e1, s1, o1, vp1, &
                          & a2, e2, s2, o2, vp2, &
                          & t0, dt_min, tf
close (20)
!------------------------------------------------------

!------------------------ START -----------------------
! Time
call cpu_time (start_time)

! Output file
open (10, file=trim (filename)//".txt", status='replace')
!------------------------------------------------------

!------------------------ LOOP ------------------------

! Initial parameters
call set_y0 (s0, o0, y0)
call set_yi (a1, e1, s1, o1, vp1, y1)
call set_yi (a2, e2, s2, o2, vp2, y2)
call set_big_y (y0, y1, y2, y)

! Aux parameters (Roche)
aux1 = (r(2) / 0.462) * (m(1) / m(2))**(1/3.)
aux2 = (r(3) / 0.462) * (m(1) / m(3))**(1/3.)

! Run variables
t       = t0         ! Init time
dt      = Logt       ! First dt
t_add   = Logt       ! For incrementing t_out
t_out   = t0 + t_add ! First output time
dt_adap = dt         ! For adaptive step
i       = 0          ! Counter

! Write initial conditions
write (*, "(I11, 3(E14.4, 1X))") i, t, dt, dt_adap
write (10,*) y, n_f (y(3), mu(2)), n_f (y(8), mu(3)), t / YR2DAY, dt / YR2DAY, dt_adap / YR2DAY

! Do loop
do while (t < tf)
    ! CHECK ENCOUNTERS
    !! s0, o0, a1, K1, s1, o1, H1, a2, K2, s2, o2, H2
    !!  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12
    if (y(3) .le. aux1) then
        write (*,*) "End of RUN. [Encounter of planet 1]"
        exit
    end if
    if (y(8) .le. aux2) then
        write (*,*) "End of RUN. [Encounter of planet 2]"
        exit
    end if
    
    !!! Execute an integration method (uncomment/edit one of theese)
    ! call implicit_caller (t, y, dt, dydtidall, euler_centred, max_iter, e_tol, dt_min, ynew)
    ! call integ_caller (t, y, dt, dydtidall, rungek6, dt_min, ynew)
    ! call rk_half_step_caller (t, y, dt_adap, dydtidall, rungek4, 4, e_tol, beta, dt_min, dt, ynew)
    call embedded_caller (t, y, dt_adap, dydtidall, Bulirsch_Stoer, e_tol, beta, dt_min, dt, ynew)
    
    !! Modulate and avoid too small angles
    ynew(2)  = max (MIN_VAL, mod (ynew(2), TWOPI))
    ynew(6)  = max (MIN_VAL, mod (ynew(6), TWOPI))
    ynew(11) = max (MIN_VAL, mod (ynew(11), TWOPI))
    
    ! Update parameters
    i  = i + 1
    t  = t + dt
    y  = ynew
    
    ! Output & dt
    if (t .ge. t_out) then
        do while (t_out .le. t)
            t_add = t_add * Logt
            t_out = t0 + t_add
        end do
        if ((t_out > tf) .or. ((abs (t_out - tf) / tf) .le. e_tol)) then
            t_out = tf
        end if
        dt = t_out - t
        write (*, "(I11, 3(E14.4, 1X))") i, t / YR2DAY, dt / YR2DAY, dt_adap / YR2DAY
        write (10,*) y, n_f (y(3), mu(2)), n_f (y(8), mu(3)), t, dt, dt_adap
    end if
    
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
aux1 = final_time - start_time
if (aux1 > 3600.) then
    print '("Total running time:",F10.4," [hs].")', aux1 / 3600.
else if (aux1 > 60.) then
    print '("Total running time:",F10.4," [min].")', aux1 / 60.
else
    print '("Total running time:",F10.4," [seg].")', aux1
end if

!------------------------------------------------------

end program tidal
