program tidal

use run
use tidall
use integrators

implicit none
!------------------ INPUT PARAMETERS ------------------

! Objects
!! Object 0
m0      = 1.0             ! [Ms]
s0      = TWOPI / 28.     ! [rad day⁻¹]
o0      = 25. * PI / 180. ! [rad]
radius0 = 695700. * KM2UA ! [UA]
alpha0  = 0.4             ! []
Q0      = 1.e6            ! [?]

!! Object 1
m1      = 1. * MJ2MS      ! [Ms]
a1      = 0.05            ! [UA] 
e1      = 0.1             ! []
s1      = TWOPI / 0.01    ! [rad day⁻¹]
o1      = 80. * PI / 180. ! [rad]
radius1 = 69911. * KM2UA  ! [UA]
alpha1  = 0.4             ! []
Q1      = 1.e5            ! [?]

! Run conditions
t0       = 0. * YR2DAY    ! [days]
dt       = 1. * YR2DAY   ! [days] ![First & min]
tf       = 1.e10 * YR2DAY ! [days]
n_points = 3500           ! N_output

! Integration conditions
beta   = 0.95 ! Learning rate
e_tol  = 1e-6 ! Approx Absolute e_calc (|Ysol - Ypred|)

! Output
filename = "Salida.txt"
!------------------------------------------------------

!-------------- DERIVED GLOBAL PARAMETERS --------------
! Calculated constants
m1p   = (m0 * m1) / (m0 + m1)     ! []
mu    = G * (m0 + m1)             ! [AU³ days⁻²]
C0    = alpha0 * m0 * radius0**2  ! [Ms AU²]
C1    = alpha1 * m1 * radius1**2  ! [Ms AU²]
K0    = Ki (a1, m1, radius0, Q0)  ! [?]
K1    = Ki (a1, m0, radius1, Q1)  ! [?]

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
write (20, '(17(A,1X))') "m0", "r0", "al0", "Q0", "m1", "r1", "al1", "Q1", &
                       & "a1", "e1", "s1", "o1", "s0", "o0", &
                       & "t0", "dt", "tf"
write (20,'(17(E16.9,1X))') m0, radius0, alpha0, Q0, &
                          & m1, radius1, alpha1, Q1, &
                          & a1, e1, s1, o1, s0, o0,  &
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
call set_y (a1, e1, s1, o1, s0, o0, y)
t       = t0
t_out   = t0 + t_add
dt_min  = dt
dt_adap = dt ! For adaptive step
i       = 0


! Do loop
do while (t < tf)
    
    ! CHECK
    if (any (((y - y) .ne. 0.d0)) .or. (y(1) < 0)) then
        write (*,*) "End of RUN. [Encounter]"
        exit
    end if
    
    ! Output    
    if (t >= t_out) then        
        t_add = t_add * Logt
        t_out = t0 + t_add
        write (*, "(I11, 8(E14.4, 1X))") i, y(1), y(2), y(3), y(4), y(5), y(6), t / YR2DAY, dt / YR2DAY
        write (10,*) y(1), y(2), y(3), y(4), y(5), y(6), t
    end if
    
    !!! Execute an integration method (uncomment one of theese)
!     call integ_caller (t, y, dt, dydtidal, rungek4, ynew)
!     call rec_rk_adap (t, y, dt_adap, dydtidal, rungek4, 4, e_tol, beta, dt_min, dt, ynew)
     call rec_rk4_5 (t, y, dt_adap, dydtidal, e_tol, beta, dt_min, dt, ynew)
    
    ! Modulate obliquity angles
    ynew(4) = mod (ynew(4), TWOPI)
    ynew(6) = mod (ynew(6), TWOPI)

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


