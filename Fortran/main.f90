program tidal

use const
use vars
use auxs
use params_f
use funcs_f
use derivates
use rk4

implicit none

!---------------- INPUT PARAMETERS ----------------

! Objects
!! Object 0
m0      = 1.0             ! [Ms]
s0      = TWOPI / 28.     ! [rad day⁻¹]
o0      = 0.0             ! [rad]
radius0 = 695700. * KM2UA ! [UA]
alpha0  = 0.4             ! []
Q0      = 1.e6            ! [?]

!! Object 1
m1      = 1. * MJ2MS      ! [Ms]
a1      = 0.05            ! [UA] 
e1      = 0.1             ! []
s1      = TWOPI / 1.      ! [rad day⁻¹]
o1      = 0.0             ! [rad]
radius1 = 69911. * KM2UA  ! [UA]
alpha1  = 0.4             ! []
Q1      = 1.e5            ! [?]

!Run conditions
t0       = 0. * YR2DAY    ! [days]
dt       = 25. * YR2DAY   ! [days]
tf       = 2.5e9 * YR2DAY ! [days]
n_points = 3000           ! N_output

!Output
filename = "Salida.txt"
!---------------------------------------------------


! Calculated constants
m1p    = (m0 * m1) / (m0 + m1)                           ! []
mu     = G * (m0 + m1)                                   ! [AU³ days⁻²]
C0     = alpha0 * m0 * radius0**2                        ! [Ms AU²]
C1     = alpha1 * m1 * radius1**2                        ! [Ms AU²]
K0cte  = 4.5 * G * m1**2 * radius0**5 / (Q0 * sqrt (mu)) ! [?]
K1cte  = 4.5 * G * m0**2 * radius1**5 / (Q1 * sqrt (mu)) ! [?]

! Calculated parameters
n1 = nf  (a1)
K0 = K0f (a1)
K1 = K1f (a1)

! Calculated LOOP parameters
n_iter = int ((tf - t0) / dt, kind=8)
Logt   = exp (log (tf - t0) / (n_points - 1))
t_add  = Logt
print *, "Approximate Iterations:", n_iter

!------------------   LOOP   -----------------------
inquire (file=filename, exist=file_exists)

if (file_exists) then
    print '("File ",A10, "already exists.")', filename
    print*, "Enter an option:"
    print*, "      R: Rewrite file. (Default)"
    print*, "      E: Exit."
    read*, selection
    select case (selection)
        case ("R")
            print*, "Rewriting..."
        case default
            print*, ("Exiting.")
            stop (0)
    end select
end if


open (10, file=filename, status='replace')
t      = t0
t_out  = t0 + t_add
call cpu_time (start_time)
do i = 0, n_iter
    t = t + dt;

    call integrk4 (a1, dadt,  a10)
    call integrk4 (e1, dedt,  e10)
    call integrk4 (s1, ds1dt, s10)
    call integrk4 (o1, do0dt, o10)
    call integrk4 (s0, ds0dt, s00)
    call integrk4 (o0, do0dt, o00)

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
            print '("Total iterations:",i11)', i
            exit
        end if
        t_add = t_add * Logt
        t_out = t0 + t_add
        ! write (*,*) a1, e1, s1, o1, s0, o0, t / YR2DAY
        write (*, "(F10.7)")  a1, e1, s1, o1, s0, o0, t / YR2DAY
        print '("Iteration:",i10)', i
        write (10,*) a1, e1, s1, o1, s0, o0, t
    end if
end do

close (10)

call cpu_time (final_time)

print '("Total running time:",f8.3," [sec].")', final_time - start_time

end program tidal

