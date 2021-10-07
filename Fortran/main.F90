module rkall

    use integrators
    use tidall

    implicit none
       
    contains
        subroutine integ_tidal (t, dt, integ, a1out, e1out, s1out, o1out, s0out, o0out)
            implicit none
            real*8, intent (in)  :: t, dt
            procedure(integ_tem) :: integ
            real*8, intent (out) :: a1out, e1out, s1out, o1out, s0out, o0out

            call integ (t, a1, dt,  dadt, a1out)
            call integ (t, e1, dt,  dedt, e1out)
            call integ (t, s1, dt, ds1dt, s1out)
            call integ (t, o1, dt, do1dt, o1out)
            call integ (t, s0, dt, ds0dt, s0out)
            call integ (t, o0, dt, do0dt, o0out)
        end subroutine integ_tidal

        subroutine rksolve_tidal (t, dt, m, a1out, e1out, s1out, o1out, s0out, o0out)
            implicit none
            real*8, intent(in)                  :: t, dt
            real*8, dimension (:,:), intent(in) :: m
            real*8, intent (out)                :: a1out, e1out, s1out, o1out, s0out, o0out

            call rksolve (t, a1, dt, dadt,  m, a1out)
            call rksolve (t, e1, dt, dedt,  m, e1out)
            call rksolve (t, s1, dt, ds1dt, m, s1out)
            call rksolve (t, o1, dt, do1dt, m, o1out)
            call rksolve (t, s0, dt, ds0dt, m, s0out)
            call rksolve (t, o0, dt, do0dt, m, o0out)
        end subroutine rksolve_tidal

        recursive subroutine rec_rk_adap_tidal (t, dt, integr, p, e_tol, beta, dt_min, a10, e10, s10, o10, s00, o00)
            implicit none
            integer, intent(in)    :: p
            procedure(integ_tem)   :: integr
            real*8, intent (in)    :: t, e_tol, beta, dt_min
            real*8, intent (inout) :: dt
            real*8, intent (out)   :: a10, e10, s10, o10, s00, o00
            real*8                 :: a11, e11, s11, o11, s01, o01
            real*8                 :: e_calc

            dt = max (dt, dt_min)
            call integ_tidal (t, dt, integr, a10, e10, s10, o10, s00, o00)
            call integ_tidal (t, dt * 0.5, integr, a11, e11, s11, o11, s01, o01)

            e_calc = norm2 ((/a10-a11, e10-e11, s10-s11, o10-o11, s00-s01, o00-o01/)) / (2.**p - 1.)
            if (e_calc < e_tol) then
                dt = max (beta * dt * (e_tol / e_calc)**(1./real (p)), dt_min)
            else
                dt = beta * dt * (e_tol / e_calc)**(1./real (p + 1))
                if ((isnan (dt)) .or. (dt <= dt_min)) then
                    dt = dt_min
                    call integ_tidal (t, dt, integr, a10, e10, s10, o10, s00, o00)
                else
                    call rec_rk_adap_tidal (t, dt, integr, p, e_tol, beta, dt_min, a10, e10, s10, o10, s00, o00)
                end if
            end if
        end subroutine rec_rk_adap_tidal

        recursive subroutine rec_rk4_5_tidal (t, dt, e_tol, beta, dt_min, a10, e10, s10, o10, s00, o00)
            implicit none
            real*8, intent (in)                :: t, e_tol, beta, dt_min
            real*8, intent (inout)             :: dt
            real*8, intent (out)               :: a10, e10, s10, o10, s00, o00
            real*8                             :: a11, e11, s11, o11, s01, o01
            real*8                             :: e_calc
            real*8, dimension(:), allocatable  :: rka, rke, rks1, rko1, rks0, rko0
            real*8, parameter, dimension(49)   :: &
            m =  (/  0.,         0.,          0.,          0.,           0.,      0.,    0., & !k1
               &   0.25,       0.25,          0.,          0.,           0.,      0.,    0., & !k2
               &  0.375,      3/32.,       9/32.,          0.,           0.,      0.,    0., & !k3
               & 12/13., 1932/2197., -7200/2197.,  7296/2197.,           0.,      0.,    0., & !k4
               &     1.,   439/216.,         -8.,   3680/513.,   -845/4104.,      0.,    0., & !k5
               &    0.5,     -8/27.,          2., -3544/2565.,   1859/4104., -11/40.,    0., & !k6
               &     0.,    16/135.,          0., 6656/12825., 28561/56430.,   -0.18, 2/55. /) !y
            real*8, parameter, dimension(6)    :: &
               & c0 = (/25/216., 0.,  1408/2565.,   2197/4104.,  -0.2,    0./), &
               & c1 = (/16/135., 0., 6656/12825., 28561/56430., -0.18, 2/55./)


            dt = max (dt, dt_min)

            call get_rks (t, a1, dt, dadt,  reshape (m, shape=(/7,7/)), rka)
            call get_rks (t, e1, dt, dedt,  reshape (m, shape=(/7,7/)), rke)
            call get_rks (t, s1, dt, ds1dt, reshape (m, shape=(/7,7/)), rks1)
            call get_rks (t, o1, dt, do1dt, reshape (m, shape=(/7,7/)), rko1)
            call get_rks (t, s0, dt, ds0dt, reshape (m, shape=(/7,7/)), rks0)
            call get_rks (t, o0, dt, do0dt, reshape (m, shape=(/7,7/)), rko0)

            a11  = a1 + dt * dot_product (c0, rka)
            e11  = e1 + dt * dot_product (c0, rke)
            s11  = s1 + dt * dot_product (c0, rks1)
            o11  = o1 + dt * dot_product (c0, rko1)
            s01  = s0 + dt * dot_product (c0, rks0)
            o01  = o0 + dt * dot_product (c0, rko0)
            a10  = a1 + dt * dot_product (c1, rka)
            e10  = e1 + dt * dot_product (c1, rke)
            s10  = s1 + dt * dot_product (c1, rks1)
            o10  = o1 + dt * dot_product (c1, rko1)
            s00  = s0 + dt * dot_product (c1, rks0)
            o00  = o0 + dt * dot_product (c1, rko0)

            e_calc = norm2 ((/a10-a11, e10-e11, s10-s11, o10-o11, s00-s01, o00-o01/))
            if (e_calc < e_tol) then
                dt = max (beta * dt * (e_tol / e_calc)**0.25, dt_min)
            else
                dt = beta * dt * (e_tol / e_calc)**0.2
                if ((isnan (dt)) .or. (dt < dt_min)) then
                    dt = dt_min
                    call rksolve_tidal (t, dt, reshape (m, (/7,7/)), a10, e10, s10, o10, s00, o00) 
                else
                    call rec_rk4_5_tidal (t, dt, e_tol, beta, dt_min, a10, e10, s10, o10, s00, o00)
                end if 
            end if
        end subroutine rec_rk4_5_tidal
       
end module rkall 


program tidal

use run
use tidall
use rkall
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
dt       = 10. * YR2DAY   ! [days] ![First & min]
tf       = 1.e10 * YR2DAY ! [days]
n_points = 3500           ! N_output

!Integration conditions
beta   = 0.95 ! Learning rate
e_tol  = 1e-7 ! Approx Absolute e_calcor (|Ysol - Ypred|)

!Output
filename = "Salida2.txt"
!---------------------------------------------------


! Calculated constants
m1p   = (m0 * m1) / (m0 + m1)     ! []
mu    = G * (m0 + m1)             ! [AU³ days⁻²]
C0    = alpha0 * m0 * radius0**2  ! [Ms AU²]
C1    = alpha1 * m1 * radius1**2  ! [Ms AU²]
K0    = Ki (a1, m1, radius0, Q0)   ! [?]
K1    = Ki (a1, m0, radius1, Q1)   ! [?]
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


! START
call cpu_time (start_time)

open (10, file=trim (filename), status='replace')
t      = t0
t_out  = t0 + t_add
dt_min = dt
i      = 0

! LOOP
do while (t < tf)

    ! Calculated parameters
    n1   = ni (a1)
    AM   = AngMom (a1, e1)
    
    call set_cosos (o0, o1, cos0, cos1)
    call set_fs (e1, fe1, fe2, fe3, fe4, fe5)
    
    ! call integ_tidal (t, dt, rungek4, a10, e10, s10, o10, s00, o00)
    call rec_rk_adap_tidal (t, dt, rungek4, 4, e_tol, beta, dt_min, a10, e10, s10, o10, s00, o00)
    ! call rec_rk4_5_tidal(t, dt, e_tol, beta, dt_min, a10, e10, s10, o10, s00, o00)

    o10 = mod (o10, TWOPI)
    o00 = mod (o00, TWOPI)

    t = t + dt
    a1 = a10
    e1 = e10
    s0 = s00
    o0 = o00
    s1 = s10    
    o1 = o10
    
    i = i + 1
    if (t >= t_out) then
        if ((isnan (a1)) .or. (a1 < 0)) then
            write (*,*) "End of RUN. [Encounter]"
            exit
        end if
        t_add = t_add * Logt
        t_out = t0 + t_add
        write (*, "(I11, 8(E14.4, 1X))") i, a1, e1, s1, o1, s0, o0, t / YR2DAY, dt / YR2DAY
        write (10,*) a1, e1, s1, o1, s0, o0, t
    end if
    
end do

print '("Total iterations:",I11)', i
print '("Final time:",E16.6," [yrs]")', t / YR2DAY

close (10)

call cpu_time (final_time)

print '("Total running time:",F10.4," [sec].")', final_time - start_time

end program tidal


