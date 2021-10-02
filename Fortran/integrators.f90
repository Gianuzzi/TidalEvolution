module integrators
    implicit none
    integer*8, parameter :: MAX_ITER = 10000
    real*8, parameter    :: MIN_ERR = 1e-9
    contains
        subroutine euler_forward (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)  :: t, y, dt
            real*8, external    :: dydt
            real*8, intent(out) :: ynew
            
            ynew = y + dydt (t, y) * dt
        end subroutine euler_forward

        subroutine euler_backward (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)  :: t, y, dt
            real*8, external    :: dydt
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
            real*8, external    :: dydt
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
            ynew = y + 0.5 * (dy1 + dydt (t, y)) * dt;
        end subroutine euler_centred

        subroutine rungek2 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)  :: t, y, dt
            real*8, external    :: dydt
            real*8, intent(out) :: ynew
            real*8              :: rk1, rk2 ! RungeKuttas
    
            rk1 = dydt (t, y)
            rk2 = dydt (t + dt, y + dt * rk1)
            ynew = y + 0.5 * (rk1 + rk2) * dt
        end subroutine rungek2

        subroutine midpoint (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)  :: t, y, dt
            real*8, external    :: dydt
            real*8, intent(out) :: ynew
            real*8              :: rk1, rk2 ! RungeKuttas
    
            rk1 = dydt (t, y)
            rk2 = dydt (t + dt * 0.5, y + dt * rk1 * 0.5)
            ynew = y + rk2 * dt
        end subroutine midpoint

        subroutine ralston (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)  :: t, y, dt
            real*8, external    :: dydt
            real*8, intent(out) :: ynew
            real*8              :: rk1, rk2 ! RungeKuttas
    
            rk1 = dydt (t, y)
            rk2 = dydt (t + dt * 0.75, y + dt * rk1 * 0.75)
            ynew = y + (rk1 + 2. * rk2) / 3. * dt
        end subroutine ralston


        subroutine rungek4 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)  :: t, y, dt
            real*8, external    :: dydt
            real*8, intent(out) :: ynew
            real*8              :: rk1, rk2, rk3, rk4 ! RungeKuttas
    
            rk1 = dydt (t, y)
            rk2 = dydt (t + 0.5 * dt, y + dt * rk1 * 0.5)
            rk3 = dydt (t + 0.5 * dt, y + dt * rk2 * 0.5)
            rk4 = dydt (t, y + dt * rk3)
            ynew = y + 1 / 6. * (rk1 + 2. * (rk2 + rk3) + rk4) * dt
        end subroutine rungek4
end module integrators