module rk4
    use const
    use derivates
    contains
        subroutine integrk4 (y, dydy, anew)
        implicit none
        real*8, intent(in)  :: y
        real*8, external    :: dydy
        real*8, intent(out) :: anew
        real*8              :: rk1, rk2, rk3, rk4 ! RungeKuttas

        rk1 = dydy (y)
        rk2 = dydy (y + dt * rk1 * 0.5)
        rk3 = dydy (y + dt * rk2 * 0.5)
        rk4 = dydy (y + dt * rk3)
        anew = y + dt * (rk1 + 2. * (rk2 + rk3) + rk4) / 6.
        end subroutine integrk4
end module rk4