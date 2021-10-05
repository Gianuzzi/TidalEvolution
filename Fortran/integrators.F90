module integrators

    implicit none
    integer*8, parameter :: MAX_ITER = 10000
    real*8, parameter    :: MIN_ERR = 1e-9

    abstract interface 

        real*8 function dydt_tem (t, y) result (der)
            real*8, intent(in)  :: t, y
        end function dydt_tem

        subroutine integ_tem (t, y, dt, dydt, ynew)
            real*8, intent(in)  :: t, y, dt
            procedure(dydt_tem) :: dydt
            real*8, intent(out) :: ynew
        end subroutine integ_tem

    end interface

    contains

        subroutine get_rks (t, y, dt, dydt, m, rk)
            implicit none
            real*8, intent(in)                                :: t, y, dt
            real*8, dimension (:,:), intent(in)               :: m
            real*8, dimension (:), allocatable, intent(inout) :: rk
            procedure(dydt_tem)                               :: dydt
            integer                                           :: N, i

            N = size (m, 1)
            allocate (rk(N-1))
            do i = 1, N - 1
                rk(i) = dydt (t + m(1,i) * dt, y + dt * dot_product (m(2:,i), rk))
            end do
        
        end subroutine get_rks

        subroutine rksolve (t, y, dt, dydt, m, ynew)
            implicit none
            real*8, intent(in)                  :: t, y, dt
            real*8, dimension (:,:), intent(in) :: m
            real*8, dimension (:), allocatable  :: rk
            procedure(dydt_tem)                 :: dydt
            real*8, intent(out)                 :: ynew
            
            call get_rks (t, y, dt, dydt, m, rk)
            ynew = y + dt * dot_product (m(2:,size (m, 1)), rk)
        end subroutine rksolve

        subroutine euler_forward (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)  :: t, y, dt
            procedure(dydt_tem) :: dydt
            real*8, intent(out) :: ynew
            
            ynew = y + dydt (t, y) * dt
        end subroutine euler_forward

        subroutine euler_backward (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)  :: t, y, dt
            procedure(dydt_tem) :: dydt
            real*8, intent(out) :: ynew
            real*8              :: y1, dy1
            integer*8           :: k

            y1 = y
            do k = 0, MAX_ITER
                dy1 = dydt (t + dt, y1) * dt
                y1  = y1 + dy1
                if (abs (dy1) <= MIN_ERR) then
                    exit
                end if 
            end do            
            ynew = y + dy1
        end subroutine euler_backward

        subroutine euler_centred (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)  :: t, y, dt
            procedure(dydt_tem) :: dydt
            real*8, intent(out) :: ynew
            real*8              :: y1, dy1
            integer*8           :: k

            y1 = y
            do k = 0, MAX_ITER
                dy1 = dydt (t + dt, y1) * dt
                y1  = y1 + dy1
                if (abs (dy1) <= MIN_ERR) then
                    exit
                end if 
            end do            
            ynew = y + 0.5 * (dydt (t + dt, y1) + dydt (t, y)) * dt;
        end subroutine euler_centred

        subroutine rungek2 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)  :: t, y, dt
            procedure(dydt_tem) :: dydt
            real*8, intent(out) :: ynew
            real*8              :: rk1, rk2 ! RungeKuttas
    
            rk1 = dydt (t, y)
            rk2 = dydt (t + dt, y + dt * rk1)
            ynew = y + 0.5 * (rk1 + rk2) * dt
        end subroutine rungek2

        subroutine midpoint2 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)  :: t, y, dt
            procedure(dydt_tem) :: dydt
            real*8, intent(out) :: ynew
            real*8              :: rk1, rk2 ! RungeKuttas
    
            rk1 = dydt (t, y)
            rk2 = dydt (t + dt * 0.5, y + dt * rk1 * 0.5)
            ynew = y + rk2 * dt
        end subroutine midpoint2

        subroutine strange2 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)  :: t, y, dt
            procedure(dydt_tem) :: dydt
            real*8, intent(out) :: ynew
            real*8              :: rk1, rk2 ! RungeKuttas
    
            rk1 = dydt (t, y)
            rk2 = dydt (t + dt * 0.75, y + dt * rk1 * 0.75)
            ynew = y + (rk1 + 2. * rk2) /3. * dt
        end subroutine strange2

        subroutine ralston2 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)  :: t, y, dt
            procedure(dydt_tem) :: dydt
            real*8, intent(out) :: ynew
            real*8              :: rk1, rk2 ! RungeKuttas
    
            rk1 = dydt (t, y)
            rk2 = dydt (t + dt * 2/3., y + dt * rk1 * 2/3.)
            ynew = y + (rk1 + 3. * rk2) * 0.25 * dt
        end subroutine ralston2

        subroutine rungek3 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)  :: t, y, dt
            procedure(dydt_tem) :: dydt
            real*8, intent(out) :: ynew
            real*8              :: rk1, rk2, rk3 ! RungeKuttas
    
            rk1 = dydt (t, y)
            rk2 = dydt (t + 0.5 * dt, y + dt * rk1 * 0.5)
            rk3 = dydt (t + dt, y + dt * (2. * rk2 - rk1))
            ynew = y + 1 /6. * (rk1 + 4. * rk2 + rk3) * dt
        end subroutine rungek3

        subroutine heun3 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)  :: t, y, dt
            procedure(dydt_tem) :: dydt
            real*8, intent(out) :: ynew
            real*8              :: rk1, rk2, rk3 ! RungeKuttas
    
            rk1 = dydt (t, y)
            rk2 = dydt (t + 1/3. * dt, y + dt * rk1 * 1/3.)
            rk3 = dydt (t + 2/3. * dt, y + dt * 2/3. * rk2)
            ynew = y + 0.25 * (rk1 + 3. * rk3) * dt
        end subroutine heun3

        subroutine ralston3 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)  :: t, y, dt
            procedure(dydt_tem) :: dydt
            real*8, intent(out) :: ynew
            real*8              :: rk1, rk2, rk3 ! RungeKuttas
    
            rk1 = dydt (t, y)
            rk2 = dydt (t + 0.5 * dt, y + dt * rk1 * 0.5)
            rk3 = dydt (t + 0.75 * dt, y + dt * 0.75 * rk2)
            ynew = y + 1/9. * (2 * rk1 + 3. * rk2 + 4. * rk3) * dt
        end subroutine ralston3

        subroutine ssprk3 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)  :: t, y, dt
            procedure(dydt_tem) :: dydt
            real*8, intent(out) :: ynew
            real*8              :: rk1, rk2, rk3 ! RungeKuttas
    
            rk1 = dydt (t, y)
            rk2 = dydt (t + dt, y + dt * rk1)
            rk3 = dydt (t + 0.5 * dt, y + dt * 0.25 * (rk1 + rk2))
            ynew = y + 1/6. * (rk1 + rk2 + 4. * rk3) * dt
        end subroutine ssprk3

        subroutine rungek4 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)  :: t, y, dt
            procedure(dydt_tem) :: dydt
            real*8, intent(out) :: ynew
            real*8              :: rk1, rk2, rk3, rk4 ! RungeKuttas
    
            rk1 = dydt (t, y)
            rk2 = dydt (t + 0.5 * dt, y + dt * rk1 * 0.5)
            rk3 = dydt (t + 0.5 * dt, y + dt * rk2 * 0.5)
            rk4 = dydt (t, y + dt * rk3)
            ynew = y + 1 /6. * (rk1 + 2. * (rk2 + rk3) + rk4) * dt
        end subroutine rungek4

        subroutine rungek4_3oct (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)  :: t, y, dt
            procedure(dydt_tem) :: dydt
            real*8, intent(out) :: ynew
            real*8              :: rk1, rk2, rk3, rk4 ! RungeKuttas

            rk1 = dydt (t, y)
            rk2 = dydt (t + 1/3. * dt, y + dt * rk1 * 1/3.)
            rk3 = dydt (t + 2/3. * dt, y + dt * (-rk2 * 1/3.))
            rk4 = dydt (t + dt, y + dt * (rk1 - rk2 + rk3))
            ynew = y + 0.125 * (rk1 + 3. * (rk2 + rk3) + rk4) * dt
        end subroutine rungek4_3oct

        subroutine Fehlberg4 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)  :: t, y, dt
            procedure(dydt_tem) :: dydt
            real*8, intent(out) :: ynew
            real*8              :: rk1, rk2, rk3, rk4, rk5, rk6 ! RungeKuttas
    
            rk1 = dydt (t, y)
            rk2 = dydt (t + 0.25 * dt, y + dt * rk1 * 0.25)
            rk3 = dydt (t + 0.375 * dt, y + dt * (3. * rk1 + 9. * rk2) /32.)
            rk4 = dydt (t + 12/13. * dt,  y + dt * (rk2 + 2. * rk3) * 0.5)
            rk5 = dydt (t + 0.75 * dt, y - dt * (0.5 * rk1 + 3. * rk4) * 0.375)
            rk6 = dydt (t, y - dt * (3. * rk1 - 2. * rk2 - 12. * (rk3 - rk4) - 8. * rk5) /7.)
            ynew = y + 1 /90. * (7. * (rk1 + rk6) + 32. * (rk2 + rk5) + 12. * rk4) * dt
        end subroutine Fehlberg4

        subroutine rungek6 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)    :: t, y, dt
            procedure(dydt_tem)   :: dydt
            real*8, intent(out)   :: ynew
            real*8, dimension(49) :: m 

            m = (/   0.,      0.,     0.,    0.,     0.,    0., 0., & !k1
                  & 0.25,   0.25,     0.,    0.,     0.,    0., 0., & !k2
                  & 0.25,  0.125,  0.125,    0.,     0.,    0., 0., & !k3
                  &  0.5,     0.,   -0.5,    1.,     0.,    0., 0., & !k4
                  & 0.75, 0.1875,     0.,    0., 0.5625,    0., 0., & !k5
                  &   1.,  -3/7.,   2/7., 12/7., -12/7.,  8/7., 0., & !k6
                  &   0.,  7/90., 16/45., 6/45., 16/45., 7/90., 0.  & !y
                  & /)
            
            call rksolve(t, y, dt, dydt, reshape(m, (/7,7/)), ynew)
        end subroutine rungek6


        recursive subroutine rec_rk4_5 (t, y, dt, dydt, e_tol, beta, dt_min, ynew)
            implicit none
            procedure(dydt_tem)                :: dydt
            real*8, intent (in)                :: y, t, e_tol, beta, dt_min
            real*8, intent (inout)             :: dt
            real*8, intent (out)               :: ynew
            real*8                             :: yaux, err
            real*8, dimension(:), allocatable  :: rk
            real*8, parameter, dimension(49)   :: m = &
            & (/  0.,         0.,         0.,          0.,           0.,      0.,    0., & !k1
            &   0.25,       0.25,         0.,          0.,           0.,      0.,    0., & !k2
            &  0.375,      3/32.,      9/32.,          0.,           0.,      0.,    0., & !k3
            & 12/13., 1932/2197., 7200/2197.,  7296/2197.,           0.,      0.,    0., & !k4
            &     1.,   439/216.,        -8.,   3680/513.,   -845/4104.,      0.,    0., & !k5
            &    0.5,     -8/27.,        -2., -3544/2565.,   1859/4104., -11/40.,    0., & !k6
            &     0.,    16/135.,         0., 6656/12825., 28561/56430.,   -0.18, 2/55. /) !y

            dt = max (dt, dt_min)
            call get_rks(t, y, dt, dydt, reshape (m, shape=(/7,7/)), rk)

            yaux = y + dt * dot_product ((/25/216., 0.,  1408/2565.,   2197/4104.,  -0.2,    0./), rk)
            ynew = y + dt * dot_product ((/16/135., 0., 6656/12825., 28561/56430., -0.18, 2/55./), rk)
            
            err = abs (ynew - yaux)
            if (err < e_tol) then
                dt = max (beta * dt * (e_tol / err)**0.25, dt_min)
            else
                dt = beta * dt * (e_tol / err)**0.2
                if ((isnan (dt)) .or. (dt < dt_min)) then
                    dt = dt_min
                    call rksolve (t, y, dt, dydt, reshape(m, (/7,7/)), ynew)
                else
                    call rec_rk4_5 (t, y, dt, dydt, e_tol, beta, dt_min, ynew)
                end if 
            end if
        end subroutine rec_rk4_5


        recursive subroutine rec_rk_adap (t, y, dt, dydt, integr, p, e_tol, beta, dt_min, ynew)
            implicit none
            integer, intent(in)    :: p
            procedure(integ_tem)   :: integr
            procedure(dydt_tem)    :: dydt
            real*8, intent (in)    :: y, t, e_tol, beta, dt_min
            real*8, intent (inout) :: dt
            real*8, intent (out)   :: ynew
            real*8                 :: yaux, err

            dt = max (dt, dt_min)
            call integr(t, y, dt,       dydt, ynew)
            call integr(t, y, dt * 0.5, dydt, yaux)

            err =  norm2 ((/ynew - yaux/)) / (2.**p - 1.)
            if (err < e_tol) then
                dt = max (beta * dt * (e_tol / err)**(1./real (p)), dt_min)
            else
                dt = beta * dt * (e_tol / err)**(1./real (p + 1))
                if ((isnan (dt)) .or. (dt <= dt_min)) then
                    dt = dt_min
                    call integr(t, y, dt, dydt, ynew)
                else
                    call rec_rk_adap (t, y, dt, dydt, integr, p, e_tol, beta, dt_min, ynew)
                end if
            end if
        end subroutine rec_rk_adap

end module integrators