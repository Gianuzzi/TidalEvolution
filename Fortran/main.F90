program tidal

use run
use tidall
use integrators

implicit none

!---------------- INPUT PARAMETERS ----------------

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

!Run conditions
t0       = 0. * YR2DAY    ! [days]
dt       = 50. * YR2DAY   ! [days]
tf       = 1.e10 * YR2DAY ! [days]
n_points = 3500           ! N_output

!Output
filename = "Salida3.txt"
!---------------------------------------------------


! Calculated constants
m1p   = (m0 * m1) / (m0 + m1)     ! []
mu    = G * (m0 + m1)             ! [AU³ days⁻²]
C0    = alpha0 * m0 * radius0**2  ! [Ms AU²]
C1    = alpha1 * m1 * radius1**2  ! [Ms AU²]
K0    = Ki(a1, m1, radius0, Q0)   ! [?]
K1    = Ki(a1, m0, radius1, Q1)   ! [?]
! Calculated LOOP parameters
n_iter = int ((tf - t0) / dt, kind=8)
Logt   = exp (log (tf - t0) / (n_points - 1))
t_add  = Logt
print *, "Approximate Iterations:", n_iter

!------------------   LOOP   -----------------------
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
                       & "a1", "e1", "s1", "s1", "s0", "o0", &
                       & "t0", "dt", "tf"
write (20,'(17(E16.9,1X))') m0, radius0, alpha0, Q0, &
                          & m1, radius1, alpha1, Q1, &
                          & a1, e1, s1, o1, s0, o0,  &
                          & t0, dt, tf
close (20)

open (10, file=trim (filename), status='replace')
t      = t0
t_out  = t0 + t_add
call cpu_time (start_time)
do i = 0, n_iter
    t = t + dt;

    ! Calculated parameters
    n1   = ni (a1)
    AM   = AngMom (a1, e1)
    
    call set_cosos(o0, o1, cos0, cos1)
    call set_fs (e1, fe1, fe2, fe3, fe4, fe5)
    
    call rungek4 (t, a1, dt, dadt,  a10)
    call rungek4 (t, e1, dt, dedt,  e10)
    call rungek4 (t, s1, dt, ds1dt, s10)
    call rungek4 (t, o1, dt, do1dt, o10)
    call rungek4 (t, s0, dt, ds0dt, s00)
    call rungek4 (t, o0, dt, do0dt, o00)

    o10 = mod (o10, TWOPI)
    o00 = mod (o00, TWOPI)

    a1 = a10
    e1 = e10
    s0 = s00
    o0 = o00
    s1 = s10    
    o1 = o10

    if (t >= t_out) then
        if ((isnan (a1)) .or. (a1 < 0)) then
            write (*,*) "End of RUN. [Encounter]"
            print '("Total iterations:",I12)', i
            exit
        end if
        t_add = t_add * Logt
        t_out = t0 + t_add
        write (*, "(I12, 7(E16.6, 1X))") i, a1, e1, s1, o1, s0, o0, t / YR2DAY
        write (10,*) a1, e1, s1, o1, s0, o0, t
    end if
end do

close (10)

call cpu_time (final_time)

print '("Total running time:",F10.4," [sec].")', final_time - start_time

end program tidal

