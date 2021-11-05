module integrators

    implicit none
    type dydt_i
        procedure(dydt_i_tem), pointer, nopass :: f_i => null ()
    end type dydt_i
    ! For adaptive step and implicit (might be overwritten)
    integer :: max_iter = 1000
    real*8  :: e_tol = 1e-10, beta = 0.9

    abstract interface

        !---------------------------------------------------------------------------------------------
        ! 1D -> 1D
        !---------------------------------------------------------------------------------------------

        ! f (t, y) = der
        real*8 function dydt_tem_1D (t, y) result (der)
            implicit none
            real*8, intent(in) :: t, y
        end function dydt_tem_1D

        ! in (t, y, dt, dydt, ynew) -> ynew
        subroutine integ_tem_1D (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)     :: t, y, dt
            real*8, intent(out)    :: ynew
            procedure(dydt_tem_1D) :: dydt
        end subroutine integ_tem_1D

        ! in (t, y, dt, dydt, ynew) -> ynew, dt
        subroutine integ_tem_adap_1D (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)     :: t, y
            real*8, intent(inout)  :: dt
            real*8, intent(out)    :: ynew
            procedure(dydt_tem_1D) :: dydt
        end subroutine integ_tem_adap_1D

        !---------------------------------------------------------------------------------------------
        ! ND -> ND
        !---------------------------------------------------------------------------------------------

        ! f_i (t, y__) = der
        real*8 function dydt_i_tem (t, y) result (der)
            implicit none
            real*8, intent(in)               :: t
            real*8, dimension(:), intent(in) :: y
        end function dydt_i_tem

        ! f__ (t, y__) = (f_i (t, y__), ..., f_n (t, y__)) = der__
        function dydt_tem (t, y) result (der)
            implicit none
            real*8, intent(in)               :: t
            real*8, dimension(:), intent(in) :: y
            real*8, dimension(size (y))      :: der
            ! Here must be every f_i defined explicitly
        end function dydt_tem

        ! in (t, y__, dt, f_i, ynew) -> ynew, dt
        subroutine integ_tem_i (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)               :: t, dt
            real*8, intent(out)              :: ynew
            real*8, dimension(:), intent(in) :: y
            procedure(dydt_i_tem)            :: dydt
        end subroutine integ_tem_i

        ! in (t, y__, dt, f__, ynew__) -> ynew__
        ! Remember that, in this case,
        !  f__ == (f_1, ..., f_N) must be
        !  pre-defined explicitly
        subroutine integ_tem (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
        end subroutine integ_tem

        ! in (t, y__, dt, f__, max_iter, e_tol, ynew__) -> ynew__
        ! Remember that, in this case,
        !  f__ == (f_1, ..., f_N) must be
        !  pre-defined explicitly
        subroutine implicit_tem (t, y, dt, dydt, max_iter, e_tol, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt, e_tol
            integer, intent(in)                      :: max_iter
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
        end subroutine implicit_tem

        ! in (t, y__, dt_adap, f__, e_tol, beta, dt_min, dt_used, ynew) -> ynew__, dt_used
        ! Remember that, in this case,
        !  f__ == (f_1, ..., f_N) must be
        !  pre-defined explicitly
        subroutine embedded_tem (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            real*8, intent(inout)                    :: dt_adap, dt_used
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
        end subroutine embedded_tem

        !---------------------------------------------------------------------------------------------
        ! WRAPPERS 
        !---------------------------------------------------------------------------------------------

        ! Here dydt_vec is a pointer to an array of (f__1, ..., f__N)
        !  so this is kind of a wrapper.
        ! dydt_tem_w (t, y__, =>f__)
        !   CREATE F (t, y__) = (f__1 (t, y__), ..., f__N (t, y__)) = der__
        !   --> (f__1 (t, y__), ..., f__N (t, y__)) 
        !   --> (der_1, ..., der_N) = der__
        function dydt_tem_w (t, y, dydt) result (der)
            import :: dydt_i
            implicit none
            real*8, intent(in)                :: t
            real*8, dimension(:), intent(in)  :: y
            type(dydt_i), dimension(size (y)) :: dydt
            real*8, dimension(size (y))       :: der
        end function dydt_tem_w

        ! Here dydt_vec is a pointer to an array of (f__1, ..., f__N);
        !  so this is kind of a wrapper.
        ! integ_tem_w (t, y__, dt, =>f__, integ, ynew__) -->
        !   CREATE F (t, y__) = (f__1 (t, y__), ..., f__N (t, y__))
        !     --> integ (t, y__, dt, F__, ynew__) --> ynew__
        subroutine integ_tem_w (t, y, dt, dydt, integ, ynew)
            import :: dydt_i
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            type(dydt_i), dimension(size (y))        :: dydt
            procedure(integ_tem)                     :: integ
            real*8, dimension(size (y)), intent(out) :: ynew
        end subroutine integ_tem_w

        ! Same as before, but for implicit_integrators
        subroutine implicit_tem_w (t, y, dt, dydt, integ, max_iter, e_tol, ynew)
            import :: dydt_i
            implicit none
            real*8, intent(in)                       :: t, dt, e_tol
            integer, intent(in)                      :: max_iter
            real*8, dimension(:), intent(in)         :: y
            type(dydt_i), dimension(size (y))        :: dydt
            procedure(implicit_tem)                  :: integ
            real*8, dimension(size (y)), intent(out) :: ynew
        end subroutine implicit_tem_w

        ! Same as before, but for embedded_integrators
        subroutine embedded_tem_w (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, integ, ynew)
            import :: dydt_i
            implicit none
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            real*8, intent(inout)                    :: dt_adap, dt_used 
            real*8, dimension(:), intent(in)         :: y
            type(dydt_i), dimension(size (y))        :: dydt
            procedure(embedded_tem)                 :: integ
            real*8, dimension(size (y)), intent(out) :: ynew
        end subroutine embedded_tem_w
        
        ! Same as before, but for rk_adap
        subroutine rk_adap_w (t, y, dt_adap, dydt, integ, p, e_tol, beta, dt_min, dt_used, ynew)
            import :: dydt_i
            integer, intent(in)                      :: p
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            procedure(integ_tem)                     :: integ
            real*8, dimension(:), intent(in)         :: y
            type(dydt_i), dimension(size (y))        :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, intent(inout)                    :: dt_adap, dt_used        
        end subroutine rk_adap_w

    end interface

    contains

        !---------------------------------------------------------------------------------------------
        ! 1D -> 1D
        !---------------------------------------------------------------------------------------------

        subroutine get_rks_1D (t, y, dt, dydt, m, rk)
            implicit none
            real*8, intent(in)                                :: t, y, dt
            real*8, dimension(:,:), intent(in)                :: m
            real*8, dimension((size (m, 1) - 1)), intent(out) :: rk
            procedure(dydt_tem_1D)                            :: dydt
            integer                                           :: i

            do i = 1, size (m, 1) - 1
                rk(i) = dydt (t + m(1,i) * dt, y + dt * dot_product (m(2:,i), rk))
            end do
        
        end subroutine get_rks_1D

        subroutine solve_rk_1D (t, y, dt, dydt, m, ynew)
            implicit none
            real*8, intent(in)                   :: t, y, dt
            real*8, dimension(:,:), intent(in)   :: m
            real*8, dimension((size (m, 1) - 1)) :: rk
            procedure(dydt_tem_1D)               :: dydt
            real*8, intent(out)                  :: ynew
            
            call get_rks_1D (t, y, dt, dydt, m, rk)
            ynew = y + dt * dot_product (m(2:,size (m, 1)), rk)
        end subroutine solve_rk_1D

        subroutine euler_forward_1D (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)     :: t, y, dt
            procedure(dydt_tem_1D) :: dydt
            real*8, intent(out)    :: ynew
            
            ynew = y + dydt (t, y) * dt
        end subroutine euler_forward_1D

        subroutine euler_backward_1D (t, y, dt, dydt, max_iter, e_tol, ynew)
            implicit none
            integer, intent(in)    :: max_iter
            real*8, intent(in)     :: t, y, e_tol, dt
            procedure(dydt_tem_1D) :: dydt
            real*8, intent(out)    :: ynew
            real*8                 :: y1, dy1
            integer*8              :: k

            y1 = y
            do k = 0, max_iter
                dy1 = dydt (t + dt, y1) * dt
                y1  = y1 + dy1
                if (abs (dy1) <= e_tol) then
                    exit
                end if 
            end do            
            ynew = y + dydt (t + dt, y1) * dt
        end subroutine euler_backward_1D

        subroutine euler_centred_1D (t, y, dt, dydt, max_iter, e_tol, ynew)
            implicit none
            integer, intent(in)    :: max_iter
            real*8, intent(in)     :: t, y, e_tol, dt
            procedure(dydt_tem_1D) :: dydt
            real*8, intent(out)    :: ynew
            real*8                 :: y1, dy1
            integer*8              :: k

            y1 = y
            do k = 0, max_iter
                dy1 = dydt (t + dt, y1) * dt
                y1  = y1 + dy1
                if (abs (dy1) <= e_tol) then
                    exit
                end if 
            end do            
            ynew = y + 0.5 * (dydt (t + dt, y1) + dydt (t, y)) * dt;
        end subroutine euler_centred_1D

        subroutine rungek2_1D (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)     :: t, y, dt
            procedure(dydt_tem_1D) :: dydt
            real*8, intent(out)    :: ynew
            real*8                 :: rk1, rk2 ! RungeKuttas

            rk1 = dydt (t, y)
            rk2 = dydt (t + dt, y + dt * rk1)
            ynew = y + 0.5 * (rk1 + rk2) * dt
        end subroutine rungek2_1D

        subroutine midpoint2_1D (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)     :: t, y, dt
            procedure(dydt_tem_1D) :: dydt
            real*8, intent(out)    :: ynew
            real*8                 :: rk1, rk2 ! RungeKuttas

            rk1 = dydt (t, y)
            rk2 = dydt (t + dt * 0.5, y + dt * rk1 * 0.5)
            ynew = y + rk2 * dt
        end subroutine midpoint2_1D

        subroutine strange2_1D (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)     :: t, y, dt
            procedure(dydt_tem_1D) :: dydt
            real*8, intent(out)    :: ynew
            real*8                 :: rk1, rk2 ! RungeKuttas

            rk1 = dydt (t, y)
            rk2 = dydt (t + dt * 0.75, y + dt * rk1 * 0.75)
            ynew = y + (rk1 + 2. * rk2) /3. * dt
        end subroutine strange2_1D

        subroutine ralston2_1D (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)     :: t, y, dt
            procedure(dydt_tem_1D) :: dydt
            real*8, intent(out)    :: ynew
            real*8                 :: rk1, rk2 ! RungeKuttas

            rk1 = dydt (t, y)
            rk2 = dydt (t + dt * 2/3., y + dt * rk1 * 2/3.)
            ynew = y + (rk1 + 3. * rk2) * 0.25 * dt
        end subroutine ralston2_1D

        subroutine rungek3_1D (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)     :: t, y, dt
            procedure(dydt_tem_1D) :: dydt
            real*8, intent(out)    :: ynew
            real*8                 :: rk1, rk2, rk3 ! RungeKuttas

            rk1 = dydt (t, y)
            rk2 = dydt (t + 0.5 * dt, y + dt * rk1 * 0.5)
            rk3 = dydt (t + dt, y + dt * (2. * rk2 - rk1))
            ynew = y + 1/6. * (rk1 + 4. * rk2 + rk3) * dt
        end subroutine rungek3_1D

        subroutine heun3_1D (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)     :: t, y, dt
            procedure(dydt_tem_1D) :: dydt
            real*8, intent(out)    :: ynew
            real*8                 :: rk1, rk2, rk3 ! RungeKuttas

            rk1 = dydt (t, y)
            rk2 = dydt (t + 1/3. * dt, y + dt * rk1 * 1/3.)
            rk3 = dydt (t + 2/3. * dt, y + dt * 2/3. * rk2)
            ynew = y + 0.25 * (rk1 + 3. * rk3) * dt
        end subroutine heun3_1D

        subroutine ralston3_1D (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)     :: t, y, dt
            procedure(dydt_tem_1D) :: dydt
            real*8, intent(out)    :: ynew
            real*8                 :: rk1, rk2, rk3 ! RungeKuttas

            rk1 = dydt (t, y)
            rk2 = dydt (t + 0.5 * dt, y + dt * rk1 * 0.5)
            rk3 = dydt (t + 0.75 * dt, y + dt * 0.75 * rk2)
            ynew = y + 1/9. * (2 * rk1 + 3. * rk2 + 4. * rk3) * dt
        end subroutine ralston3_1D

        subroutine ssprk3_1D (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)     :: t, y, dt
            procedure(dydt_tem_1D) :: dydt
            real*8, intent(out)    :: ynew
            real*8                 :: rk1, rk2, rk3 ! RungeKuttas

            rk1 = dydt (t, y)
            rk2 = dydt (t + dt, y + dt * rk1)
            rk3 = dydt (t + 0.5 * dt, y + dt * 0.25 * (rk1 + rk2))
            ynew = y + 1/6. * (rk1 + rk2 + 4. * rk3) * dt
        end subroutine ssprk3_1D

        subroutine rungek4_1D (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)     :: t, y, dt
            procedure(dydt_tem_1D) :: dydt
            real*8, intent(out)    :: ynew
            real*8                 :: rk1, rk2, rk3, rk4 ! RungeKuttas

            rk1 = dydt (t, y)
            rk2 = dydt (t + 0.5 * dt, y + dt * rk1 * 0.5)
            rk3 = dydt (t + 0.5 * dt, y + dt * rk2 * 0.5)
            rk4 = dydt (t, y + dt * rk3)
            ynew = y + 1/6. * (rk1 + 2. * (rk2 + rk3) + rk4) * dt
        end subroutine rungek4_1D

        subroutine rungek4_3oct_1D (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)     :: t, y, dt
            procedure(dydt_tem_1D) :: dydt
            real*8, intent(out)    :: ynew
            real*8                 :: rk1, rk2, rk3, rk4 ! RungeKuttas

            rk1 = dydt (t, y)
            rk2 = dydt (t + 1/3. * dt, y + dt * rk1 * 1/3.)
            rk3 = dydt (t + 2/3. * dt, y + dt * (-rk2 * 1/3.))
            rk4 = dydt (t + dt, y + dt * (rk1 - rk2 + rk3))
            ynew = y + 0.125 * (rk1 + 3. * (rk2 + rk3) + rk4) * dt
        end subroutine rungek4_3oct_1D

        subroutine Fehlberg4_1D (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)     :: t, y, dt
            procedure(dydt_tem_1D) :: dydt
            real*8, intent(out)    :: ynew
            real*8                 :: rk1, rk2, rk3, rk4, rk5, rk6 ! RungeKuttas

            rk1 = dydt (t, y)
            rk2 = dydt (t + 0.25 * dt, y + dt * rk1 * 0.25)
            rk3 = dydt (t + 0.375 * dt, y + dt * (3. * rk1 + 9. * rk2) /32.)
            rk4 = dydt (t + 12/13. * dt,  y + dt * (rk2 + 2. * rk3) * 0.5)
            rk5 = dydt (t + 0.75 * dt, y - dt * (0.5 * rk1 + 3. * rk4) * 0.375)
            rk6 = dydt (t, y - dt * (3. * rk1 - 2. * rk2 - 12. * (rk3 - rk4) - 8. * rk5) /7.)
            ynew = y + 1/90. * (7. * (rk1 + rk6) + 32. * (rk2 + rk5) + 12. * rk4) * dt
        end subroutine Fehlberg4_1D

        subroutine rungek6_1D (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)     :: t, y, dt
            procedure(dydt_tem_1D) :: dydt
            real*8, intent(out)    :: ynew
            real*8, dimension(49)  :: m 

            m = (/   0.,      0.,     0.,    0.,     0.,    0., 0., & !k1
                & 0.25,   0.25,     0.,    0.,     0.,    0., 0., & !k2
                & 0.25,  0.125,  0.125,    0.,     0.,    0., 0., & !k3
                &  0.5,     0.,   -0.5,    1.,     0.,    0., 0., & !k4
                & 0.75, 0.1875,     0.,    0., 0.5625,    0., 0., & !k5
                &   1.,  -3/7.,   2/7., 12/7., -12/7.,  8/7., 0., & !k6
                &   0.,  7/90., 16/45., 6/45., 16/45., 7/90., 0.  & !y
                & /)
            
            call solve_rk_1D (t, y, dt, dydt, reshape (m, (/7,7/)), ynew)
        end subroutine rungek6_1D

        recursive subroutine Fehlberg4_5_1D (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            procedure(dydt_tem_1D)           :: dydt
            real*8, intent(in)               :: y, t, e_tol, beta, dt_min
            real*8, intent(inout)            :: dt_adap, dt_used
            real*8, intent(out)              :: ynew
            real*8                           :: yaux, e_calc
            real*8, dimension(6)             :: rk
            real*8, parameter, dimension(49) :: m = &
            & (/  0.,         0.,          0.,          0.,           0.,      0.,    0., & !k1
            &   0.25,       0.25,          0.,          0.,           0.,      0.,    0., & !k2
            &  0.375,      3/32.,       9/32.,          0.,           0.,      0.,    0., & !k3
            & 12/13., 1932/2197., -7200/2197.,  7296/2197.,           0.,      0.,    0., & !k4
            &     1.,   439/216.,         -8.,   3680/513.,   -845/4104.,      0.,    0., & !k5
            &    0.5,     -8/27.,          2., -3544/2565.,   1859/4104., -11/40.,    0., & !k6
            &     0.,    16/135.,          0., 6656/12825., 28561/56430.,   -0.18, 2/55. /) !y

            dt_adap = max (dt_adap, dt_min)
            call get_rks_1D (t, y, dt_adap, dydt, reshape (m, shape=(/7,7/)), rk)

            yaux = y + dt_adap * dot_product ((/25/216., 0.,  1408/2565.,   2197/4104.,  -0.2,    0./), rk)
            ynew = y + dt_adap * dot_product ((/16/135., 0., 6656/12825., 28561/56430., -0.18, 2/55./), rk)
            
            e_calc = abs (ynew - yaux)
            if (e_calc < e_tol) then
                dt_used = dt_adap
                dt_adap = max (min (beta * dt_adap * (e_tol / e_calc)**0.25, dt_adap * 2.), dt_min)
            else
                dt_adap = beta * dt_adap * (e_tol / e_calc)**0.2
                if ((isnan (dt_adap)) .or. (dt_adap < dt_min)) then
                    dt_used = dt_min
                    dt_adap = dt_min
                    call solve_rk_1D (t, y, dt_adap, dydt, reshape (m, (/7,7/)), ynew)
                else
                    call Fehlberg4_5_1D (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, ynew)
                end if 
            end if
        end subroutine Fehlberg4_5_1D

        recursive subroutine rk_adap_1D (t, y, dt_adap, dydt, integ, p, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            integer, intent(in)     :: p
            procedure(integ_tem_1D) :: integ
            procedure(dydt_tem_1D)  :: dydt
            real*8, intent(in)      :: t, y, e_tol, beta, dt_min
            real*8, intent(inout)   :: dt_adap, dt_used
            real*8, intent(out)     :: ynew
            real*8                  :: yaux, e_calc

            dt_adap  = max (dt_adap, dt_min)
            call integ (t, y,       dt_adap, dydt, ynew)
            call integ (t, y, 0.5 * dt_adap, dydt, yaux)

            e_calc =  norm2 ((/ynew - yaux/)) / (2.**p - 1.)
            if (e_calc < e_tol) then
                dt_used = dt_adap
                dt_adap = max (min (beta * dt_adap * (e_tol / e_calc)**(1./real (p)), dt_adap * 2.), dt_min)
            else
                dt_adap = beta * dt_adap * (e_tol / e_calc)**(1./real (p + 1))
                if ((isnan (dt_adap)) .or. (dt_adap <= dt_min)) then
                    dt_used = dt_min
                    dt_adap = dt_min
                    call integ (t, y, dt_adap, dydt, ynew)
                else
                    call rk_adap_1D (t, y, dt_adap, dydt, integ, p, e_tol, beta, dt_min, dt_used, ynew)
                end if
            end if
        end subroutine rk_adap_1D

        !---------------------------------------------------------------------------------------------
        ! ND -> ND
        !---------------------------------------------------------------------------------------------

        ! IMPLICIT

        !! Solver
        subroutine solve_implicit (t, y, dt, dydt, max_iter, e_tol, y1)
            implicit none
            integer, intent(in)                      :: max_iter
            real*8, intent(in)                       :: t, dt, e_tol
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: y1
            real*8, dimension(size (y))              :: dy1
            integer                                  :: i

            y1 = y
            do i = 1, max_iter
                dy1 = dydt (t + dt, y1) * dt
                y1  = y1  + dy1
                if (norm2 (dy1) <= e_tol) then
                    exit
                end if
            end do
        end subroutine solve_implicit

        subroutine euler_backward (t, y, dt, dydt, max_iter, e_tol, ynew)
            implicit none
            integer, intent(in)                      :: max_iter
            real*8, intent(in)                       :: t, dt, e_tol
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: y1
            
            call solve_implicit (t, y, dt, dydt, max_iter, e_tol, y1)

            ynew = y + 0.5 * dydt (t + dt, y1) * dt
        end subroutine euler_backward

        subroutine euler_centred (t, y, dt, dydt, max_iter, e_tol, ynew)
            implicit none
            integer, intent(in)                      :: max_iter
            real*8, intent(in)                       :: t, dt, e_tol
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: y1
            
            call solve_implicit (t, y, dt, dydt, max_iter, e_tol, y1)

            ynew = y + 0.5 * (dydt (t + dt, y1) + dydt (t, y)) * dt
        end subroutine euler_centred

        ! RKs

        !! Solver
        subroutine get_rks (t, y, dt, dydt, m, rk)
            implicit none
            real*8, intent(in)                                         :: t, dt
            real*8, dimension(:), intent(in)                           :: y
            real*8, dimension(:,:), intent(in)                         :: m
            real*8, dimension((size (m,1) - 1), size (y)), intent(out) :: rk ! In columns
            real*8, dimension(size (y))                                :: rkaux
            procedure(dydt_tem)                                        :: dydt
            integer                                                    :: i, j
            
            do i = 1, size (m, 1) - 1 ! Rows
                rkaux = 0.
                do j = 1, i ! Cols (<i bc its inf triang)    
                    rkaux = rkaux + m(1+j,i) * rk(j,:)                    
                end do
                rk(i, :) = dydt (t + m(1,i) * dt, y + dt * rkaux)
            end do
        end subroutine get_rks

        !! Solver
        subroutine solve_rk (t, y, dt, dydt, m, ynew)
            implicit none
            real*8, intent(in)                            :: t, dt
            real*8, dimension(:), intent(in)              :: y
            real*8, dimension(:,:), intent(in)            :: m
            real*8, dimension((size (m,1) - 1), size (y)) :: rk
            procedure(dydt_tem)                           :: dydt
            real*8, dimension(size (y)), intent(out)      :: ynew
            integer                                       :: i 
            
            call get_rks (t, y, dt, dydt, m, rk)
            
            do i = 1, size (y)
                ynew(i) = y(i) + dt * dot_product ((/m(2:, size (m, 1))/), rk(:,i))
            end do
        end subroutine solve_rk        

        subroutine euler_forward (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(4), parameter          :: m = &
              & (/0., 0., & !k1
              &   0., 1.  & !y
              & /)
            
            call solve_rk (t, y, dt, dydt, reshape (m, (/2,2/)), ynew)
        end subroutine euler_forward

        subroutine rungek2 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(9), parameter          :: m = &
              & (/0.,  0.,  0., & !k1
              &   1.,  1.,  0., & !k2
              &   0., 0.5, 0.5  & !y
              & /)
            
            call solve_rk (t, y, dt, dydt, reshape (m, (/3,3/)), ynew)
        end subroutine rungek2

        subroutine midpoint2 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(9), parameter          :: m = &
              & (/ 0.,  0., 0., & !k1
              &   0.5, 0.5, 0., & !k2
              &    0.,  0., 1.  & !y
              & /)
            
            call solve_rk (t, y, dt, dydt, reshape (m, (/3,3/)), ynew)
        end subroutine midpoint2

        subroutine strange2 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(9), parameter          :: m = &

              & (/  0.,   0.,   0., & !k1
              &   0.75, 0.75,   0., & !k2
              &     0., 1/3., 2/3.  & !y
              & /)
            
            call solve_rk (t, y, dt, dydt, reshape (m, (/3,3/)), ynew)
        end subroutine strange2

        subroutine ralston2 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(9), parameter          :: m = &

              & (/  0.,   0.,   0., & !k1
              &   2/3., 2/3.,   0., & !k2
              &     0., 0.25, 0.75  & !y
              & /)
            
            call solve_rk (t, y, dt, dydt, reshape (m, (/3,3/)), ynew)
        end subroutine ralston2

        subroutine rungek3 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(16), parameter         :: m = &

              & (/ 0.,   0.,   0.,   0., & !k1
              &   0.5,  0.5,   0.,   0., & !k2
              &    1.,  -1.,   2.,   0., & !k3
              &    0., 1/6., 2/3., 1/6.  & !y
              & /)
            
            call solve_rk (t, y, dt, dydt, reshape (m, (/4,4/)), ynew)
        end subroutine rungek3

        subroutine heun3 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(16), parameter         :: m = &

              & (/  0.,   0.,   0.,   0., & !k1
              &   1/3., 1/3.,   0.,   0., & !k2
              &   2/3.,   0., 1/3.,   0., & !k3
              &     0., 0.25, 1/3., 0.75  & !y
              & /)
            
            call solve_rk (t, y, dt, dydt, reshape (m, (/4,4/)), ynew)
        end subroutine heun3

        subroutine ralston3 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(16), parameter         :: m = &

              & (/  0.,   0.,   0.,   0., & !k1
              &    0.5,  0.5,   0.,   0., & !k2
              &   0.75,   0., 0.75,   0., & !k3
              &     0., 2/9., 1/3., 4/9.  & !y
              & /)
            
            call solve_rk (t, y, dt, dydt, reshape (m, (/4,4/)), ynew)
        end subroutine ralston3

        subroutine ssprk3 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(16), parameter         :: m = &

              & (/ 0.,   0.,   0.,   0., & !k1
              &    1.,   1.,   0.,   0., & !k2
              &   0.5, 0.25, 0.25,   0., & !k3
              &    0., 1/6., 1/6., 2/3.  & !y
              & /)
            
            call solve_rk (t, y, dt, dydt, reshape (m, (/4,4/)), ynew)
        end subroutine ssprk3

        subroutine ralston4 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(25), parameter         :: m = &

              & (/        0.,         0.,          0.,         0.,          0., & !k1
              &          0.4,        0.4,          0.,         0.,          0., & !k2
              &   0.45573725, 0.29697761,  0.15875964,         0.,          0., & !k3
              &           1., 0.21810040, -3.05096516, 3.83286476,          0., & !k4
              &           0., 0.17476028, -0.55148066, 1.20553560, 0.017118478  & !y
              & /)
            
            call solve_rk (t, y, dt, dydt, reshape (m, (/5,5/)), ynew)
        end subroutine ralston4

        subroutine rungek4 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(25), parameter         :: m = &
              & (/ 0.,   0.,   0.,   0.,   0., & !k1
              &   0.5,  0.5,   0.,   0.,   0., & !k2
              &   0.5,   0.,  0.5,   0.,   0., & !k3
              &    1.,   0.,   0.,   1.,   0., & !k4
              &    0., 1/6., 1/3., 1/3., 1/6.  & !y
              & /)
            
            call solve_rk (t, y, dt, dydt, reshape (m, (/5,5/)), ynew)
        end subroutine rungek4

        subroutine rungek4_3oct (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(25), parameter         :: m = &

              & (/  0.,    0.,    0.,    0.,    0., & !k1
              &   1/3.,  1/3.,    0.,    0.,    0., & !k2
              &   2/3., -1/3.,    1.,    0.,    0., & !k3
              &     1.,    1.,   -1.,    1.,    0., & !k4
              &     0., 0.125, 0.375, 0.375, 0.125  & !y
              & /)
            
            call solve_rk (t, y, dt, dydt, reshape (m, (/5,5/)), ynew)
        end subroutine rungek4_3oct

        subroutine rungek6 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(49), parameter         :: m = &

               & (/   0.,     0.,     0.,    0.,     0.,    0., 0., & !k1
               & 0.25,   0.25,     0.,    0.,     0.,    0., 0., & !k2
               & 0.25,  0.125,  0.125,    0.,     0.,    0., 0., & !k3
               &  0.5,     0.,   -0.5,    1.,     0.,    0., 0., & !k4
               & 0.75, 0.1875,     0.,    0., 0.5625,    0., 0., & !k5
               &   1.,  -3/7.,   2/7., 12/7., -12/7.,  8/7., 0., & !k6
               &   0.,  7/90., 16/45., 6/45., 16/45., 7/90., 0.  & !y
               & /)
            
            call solve_rk (t, y, dt, dydt, reshape (m, (/7,7/)), ynew)
        end subroutine rungek6

        ! EMBEDDED

        !! Solver
        recursive subroutine solve_embeed (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, m, maux, o1, o2, ynew)
            implicit none
            real*8, intent(in)                            :: t, e_tol, beta, dt_min
            real*8, intent(inout)                         :: dt_adap, dt_used
            real*8, dimension(:), intent(in)              :: y
            real*8, dimension(:,:), intent(in)            :: m
            integer*4, intent(in)                         :: o1, o2
            integer*4                                     :: i, s
            real*8, dimension(size (m,1) - 1), intent(in) :: maux
            real*8, dimension(size (y))                   :: yaux
            real*8, dimension(size (m,1) - 1, size (y))   :: rk
            procedure(dydt_tem)                           :: dydt
            real*8, dimension(size (y)), intent(out)      :: ynew            
            real*8                                        :: e_calc

            s = size (m,1)

            dt_adap = max (dt_adap, dt_min)

            call get_rks (t, y, dt_adap, dydt, m, rk)
            
            do i = 1, size (y)
                yaux(i) = y(i) + dt_adap * dot_product (   maux, rk(:,i))
                ynew(i) = y(i) + dt_adap * dot_product (m(2:,s), rk(:,i))
            end do

            e_calc = norm2 (ynew - yaux)
            if (e_calc < e_tol) then
                dt_used = dt_adap
                dt_adap = max (dt_adap * min (beta * (e_tol / e_calc)**(1./o1), 2.), dt_min)
            else
                dt_adap = beta * dt_adap * (e_tol / e_calc)**(1./o2)
                if ((isnan (dt_adap)) .or. (dt_adap < dt_min)) then
                    dt_adap = dt_min
                    dt_used = dt_min
                    call solve_rk (t, y, dt_adap, dydt, m, ynew)
                else
                    call solve_embeed (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, m, maux, o1, o2, ynew)
                end if 
            end if
        end subroutine solve_embeed

        subroutine Heun_Euler1_2 (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            real*8, intent(inout)                    :: dt_adap, dt_used
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, parameter, dimension(2)          :: maux = (/ 1., 0./)
            real*8, parameter, dimension(9)          :: m    = &
               & (/0.,  0.,  0., & !k1
               &   1.,  1.,  0., & !k2
               &   0., 0.5, 0.5  & !y
               & /)

            call solve_embeed (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, reshape (m, shape=(/3,3/)), maux, 1, 2, ynew)
        end subroutine Heun_Euler1_2

        subroutine Fehlberg1_2 (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            real*8, intent(inout)                    :: dt_adap, dt_used
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, parameter, dimension(3)          :: maux = (/1/256., 255/256., 0./)
            real*8, parameter, dimension(16)         :: m    = &
               & (/0.,     0.,       0.,     0., & !k1
               &  0.5,    0.5,       0.,     0., & !k2
               &   1., 1/256., 255/256.,     0., & !k3
               &   0., 1/512., 255/256., 1/512.  & !y
               & /)

            call solve_embeed (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, reshape (m, shape=(/4,4/)), maux, 1, 2, ynew)
        end subroutine Fehlberg1_2

        subroutine Bogacki_Shampine2_3 (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            real*8, intent(inout)                    :: dt_adap, dt_used
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, parameter, dimension(4)          :: maux = (/7/24., 0.25, 1/3., 0.125/)
            real*8, parameter, dimension(25)         :: m    = &
               & (/ 0.,   0.,   0.,   0., 0., & !k1
               &   0.5,  0.5,   0.,   0., 0., & !k2
               &  0.75,   0., 0.75,   0., 0., & !k3
               &    1., 2/9., 1/3., 4/9., 0., & !k4
               &    0., 2/9., 1/3., 4/9., 0.  & !y
               & /)
               
            call solve_embeed (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, reshape (m, shape=(/5,5/)), maux, 2, 3, ynew)
        end subroutine Bogacki_Shampine2_3

        subroutine Zonneveld3_4 (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            real*8, intent(inout)                    :: dt_adap, dt_used
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, parameter, dimension(5)          :: maux = (/-0.5, 7/3., 7/3., 13/16., -16/3./)
            real*8, parameter, dimension(36)         :: m    = &
               & (/ 0.,    0.,    0.,     0.,     0.,  0., & !k1
               &   0.5,   0.5,    0.,     0.,     0.,  0., & !k2
               &   0.5,    0.,   0.5,     0.,     0.,  0., & !k3
               &    1.,    0.,    0.,     1.,     0.,  0., & !k4
               &  0.75, 5/32., 7/32., 13/32., -1/32.,  0., & !k5
               &    0.,  1/6.,  1/3.,   1/3.,   1/6.,  0.  & !y
               & /)
            call solve_embeed (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, reshape (m, shape=(/6,6/)), maux, 3, 4, ynew)
        end subroutine Zonneveld3_4

        subroutine Merson4_5 (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            real*8, intent(inout)                    :: dt_adap, dt_used
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, parameter, dimension(5)          :: maux = (/1/6., 0., 0., 2/3., 1/6./)
            real*8, parameter, dimension(36)         :: m    = &
               & (/ 0.,    0.,   0.,     0.,  0.,  0., & !k1
               &  1/3.,  1/3.,   0.,     0.,  0.,  0., & !k2
               &  1/3.,  1/6., 1/6.,     0.,  0.,  0., & !k3
               &   0.5, 0.125,   0.,  0.375,  0.,  0., & !k4
               &    1.,   0.5,   0.,   -1.5,  2.,  0., & !k5
               &    0.,   0.1,   0.,    0.3, 0.4, 0.2  & !y
               & /)
            call solve_embeed (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, reshape (m, shape=(/6,6/)), maux, 4, 5, ynew)
        end subroutine Merson4_5

        subroutine Fehlberg4_5 (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            real*8, intent(inout)                    :: dt_adap, dt_used
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, parameter, dimension(6)          :: maux = (/25/216., 0., 1408/2565., 2197/4104., -0.2, 0./)
            real*8, parameter, dimension(49)         :: m    = &
               & (/  0.,         0.,          0.,          0.,           0.,      0.,    0., & !k1
               &   0.25,       0.25,          0.,          0.,           0.,      0.,    0., & !k2
               &  0.375,      3/32.,       9/32.,          0.,           0.,      0.,    0., & !k3
               & 12/13., 1932/2197., -7200/2197.,  7296/2197.,           0.,      0.,    0., & !k4
               &     1.,   439/216.,         -8.,   3680/513.,   -845/4104.,      0.,    0., & !k5
               &    0.5,     -8/27.,          2., -3544/2565.,   1859/4104., -11/40.,    0., & !k6
               &     0.,    16/135.,          0., 6656/12825., 28561/56430.,   -0.18, 2/55.  & !y
               & /)

            call solve_embeed (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, reshape (m, shape=(/7,7/)), maux, 4, 5, ynew)
        end subroutine Fehlberg4_5

        subroutine Cash_Karp4_5 (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            real*8, intent(inout)                    :: dt_adap, dt_used
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, parameter, dimension(6)          :: maux = (/2825/27648., 0., 18575/48384., 13525/55296., 277/14336., 0.25/)
            real*8, parameter, dimension(49)         :: m    = &
               & (/  0.,          0.,       0.,         0.,            0.,        0.,        0., & !k1
               &    0.2,         0.2,       0.,         0.,            0.,        0.,        0., & !k2
               &    0.3,       0.075,    0.225,         0.,            0.,        0.,        0., & !k3
               &    0.6,         0.3,     -0.9,        1.2,            0.,        0.,        0., & !k4
               &     1.,     -11/54.,      2.5,    -70/27.,        35/27.,        0.,        0., & !k5
               &  0.875, 1631/55296., 175/512., 575/13824., 44275/110592., 253/4096.,        0., & !k6
               &     0.,     37/378.,       0.,   250/621.,      125/594.,        0., 512/1771.  & !y
               & /)

            call solve_embeed (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, reshape (m, shape=(/7,7/)), maux, 4, 5, ynew)
        end subroutine Cash_Karp4_5

        subroutine Dormand_Prince4_5 (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            real*8, intent(inout)                    :: dt_adap, dt_used
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, parameter, dimension(7)          :: maux = &
               & (/5179/57600., 0., 7571/16695., 393/640., -92097/339200., 187/2100., 0.025/)
            real*8, parameter, dimension(64)         :: m    = &
               & (/  0.,          0.,           0.,          0.,        0.,           0.,     0., 0., & !k1
               &    0.2,         0.2,           0.,          0.,        0.,           0.,     0., 0., & !k2
               &    0.3,       0.075,        0.225,          0.,        0.,           0.,     0., 0., & !k3
               &    0.8,      44/45.,      -56/15.,       32/9.,        0.,           0.,     0., 0., & !k4
               &   8/9., 19372/6561., -25360/2187., 64448/6561., -212/729.,           0.,     0., 0., & !k5
               &     1.,  9017/3168.,     -355/33., 46732/5247.,   49/176., -5103/18656.,     0., 0., & !k6
               &     1.,     35/384.,           0.,   500/1113.,  125/192.,   2187/6784., 11/84., 0., & !k7
               &     0.,     35/384.,           0.,   500/1113.,  125/192.,  -2187/6784., 11/84., 0.  & !y
               & /)

            call solve_embeed (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, reshape (m, shape=(/8,8/)), maux, 4, 5, ynew)
        end subroutine Dormand_Prince4_5

        subroutine Verner5_6 (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            real*8, intent(inout)                    :: dt_adap, dt_used
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, parameter, dimension(8)          :: maux = (/13/160., 0., 2375/5984., 5/16., 12/85., 3/44., 0., 0./)
            real*8, parameter, dimension(81)         :: m    = &
               & (/  0.,           0.,       0.,            0.,         0.,           0., 0.,          0.,      0., & !k1
               &   1/6.,         1/6.,       0.,            0.,         0.,           0., 0.,          0.,      0., & !k2
               &  4/15.,        4/75.,   16/75.,            0.,         0.,           0., 0.,          0.,      0., & !k3
               &   2/3.,         5/6.,    -8/3.,           2.5,         0.,           0., 0.,          0.,      0., & !k4
               &   5/6.,     -165/64.,    55/6.,      -425/64.,     85/96.,           0., 0.,          0.,      0., & !k5
               &     1.,          2.4,      -8.,     4015/612.,    -11/36.,      88/255., 0.,          0.,      0., & !k6
               &  1/15., -8263/15000.,  124/75.,     -643/680.,     -0.324,  2484/10625., 0.,          0.,      0., & !k8
               &     1.,   3501/1720., -300/43., 297275/52632., -319/2322., 24068/84065., 0., 3850/26703.,      0., & !k8
               &     0.,        0.075,       0.,     875/2244.,     23/72.,    264/1955., 0.,  125/11592., 43/616.  & !y
               & /)
         
            call solve_embeed (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, reshape (m, shape=(/9,9/)), maux, 5, 6, ynew)
        end subroutine Verner5_6

        subroutine Fehlberg7_8 (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            real*8, intent(inout)                    :: dt_adap, dt_used
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, parameter, dimension(13)         :: maux = &
               & (/41/840., 0., 0., 0., 0., 34/105., 9/35., 9/35., 9/280., 9/280., 41/840., 0., 0./)
            real*8, parameter, dimension(196)        :: m    = &
               & (/  0.,          0.,    0.,      0.,        0.,         0.,       0.,         0.,     0.,      0.,     0., &
                         & 0.,      0.,      0., & !k1
               &  2/27.,       2/27.,    0.,      0.,        0.,         0.,       0.,         0.,     0.,      0.,     0., &
                         & 0.,      0.,      0., & !k2
               &   1/9.,       1/36., 1/12.,      0.,        0.,         0.,       0.,         0.,     0.,      0.,     0., &
                         & 0.,      0.,      0., & !k3
               &   1/6.,       1/24.,    0.,   0.125,        0.,         0.,       0.,         0.,     0.,      0.,     0., &
                         & 0.,      0.,      0., & !k4
               &  5/12.,       5/12.,    0., -1.5625,    1.5625,         0.,       0.,         0.,     0.,      0.,     0., &
                         & 0.,      0.,      0., & !k5
               &    0.5,        0.05,    0.,      0.,      0.25,        0.2,       0.,         0.,     0.,      0.,     0., &
                         & 0.,      0.,      0., & !k6
               &   5/6.,    -25/108.,    0.,      0.,  125/108.,    -65/27.,  125/54.,         0.,     0.,      0.,     0., &
                         & 0.,      0.,      0., & !k7
               &   1/6.,     31/300.,    0.,      0.,        0.,    61/225.,    -2/9.,    13/900.,     0.,      0.,     0., &
                         & 0.,      0.,      0., & !k8
               &   2/3.,          2.,    0.,      0.,    -53/6.,    704/45.,  -107/9.,     67/90.,     3.,      0.,     0., &
                         & 0.,      0.,      0., & !k9
               &   1/3.,    -91/108.,    0.,      0.,   23/108.,  -976/135.,  311/54.,    -19/60.,  17/6.,  -1/12.,     0., &
                         & 0.,      0.,      0., & !k10
               &     1.,  2383/4100.,    0.,      0., -341/164., 4496/1025., -301/82., 2133/4100., 45/82., 45/164., 18/41., &
                         & 0.,      0.,      0., & !k11
               &     0.,      3/205.,    0.,      0.,        0.,         0.,   -6/41.,    -3/205., -3/41.,   3/41.,  6/41., &
                         & 0.,      0.,      0., & !k12
               &     1., -1777/4100.,    0.,      0., -341/164., 4496/1025., -289/82., 2193/4100., 51/82., 33/164., 19/41., &
                         & 0.,      1.,      0., & !k13
               &     0.,          0.,    0.,      0.,        0.,         0.,  34/105.,      9/35.,  9/35.,  9/280., 9/280., &
                         & 0., 41/840., 41/840.  & !y
               & /)

            call solve_embeed (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, reshape (m, shape=(/14,14/)), maux, 7, 8, ynew)
        end subroutine Fehlberg7_8

        subroutine Dormand_Prince7_8 (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            real*8, intent(inout)                    :: dt_adap, dt_used
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, parameter, dimension(13)         :: maux = &
               & (/14005451./335480064., 0., 0., 0., 0., -59238493./1068277825., 181606767./758867731., 561292985./797845732.,&
                         & -1041891430./1371343529., 760417239./1151165299., 118820643./751138087., -528747749./2220607170., 0.25/)
            real*8, parameter, dimension(196)        :: m    = &
               & (/                   0.,                      0.,     0.,        0.,                         0., &
                         &                      0.,                        0.,                        0., &
                         &                        0.,                        0.,                      0., & 
                         &                     0.,    0., 0., & !k1
               &                   1/18.,                   1/18.,     0.,        0.,                         0., &
                         &                      0.,                        0.,                        0., &
                         &                        0.,                        0.,                      0., & 
                         &                     0.,    0., 0., & !k2
               &                   1/12.,                   1/48., 0.0625,        0.,                         0., &
                         &                      0.,                        0.,                        0., &
                         &                        0.,                        0.,                      0., & 
                         &                     0.,    0., 0., & !k3
               &                   0.125,                 0.03125,     0.,   0.09375,                         0., &
                         &                      0.,                        0.,                        0., &
                         &                        0.,                        0.,                      0., & 
                         &                     0.,    0., 0., & !k4
               &                  0.3125,                  0.3125,     0., -1.171875,                   1.171875, &
                         &                      0.,                        0.,                        0., &
                         &                        0.,                        0.,                      0., & 
                         &                     0.,    0., 0., & !k5
               &                   0.375,                  0.0375,     0.,        0.,                     0.1875, &
                         &                    0.15,                        0.,                        0., &
                         &                        0.,                        0.,                      0., & 
                         &                     0.,    0., 0., & !k6
               &                  0.1475,    29443841./614563906.,     0.,        0.,       77736538./692538347., &
                         &  -28693883./1125000000.,     23124283./1800000000.,                        0., &
                         &                        0.,                        0.,                      0., & 
                         &                     0.,    0., 0., & !k7
               &                   0.465,    16016141./946692911.,     0.,        0.,       61564180./158732637., &
                         &    22789713./633445777.,    545815736./2771057229.,   -180193667./1043307555., &
                         &                        0.,                        0.,                      0., & 
                         &                     0.,    0., 0., & !k8
               & 5490023248./9719169821.,    39632708./573591083.,     0.,        0.,     -433636366./683701615., &
                         & -421739975./2616292301.,     100302831./723423059.,     790204164./839813087., &
                         &    800635310./3783071287.,                        0.,                      0., & 
                         &                     0.,    0., 0., & !k9
               &                    0.65,  246121993./1340847787.,     0.,        0., -37695042795./15268766246., &
                         & -309121744./1061227803.,     -12992083./490766935.,   6005943493./2108947869., &
                         &    393006217./1396673457.,    123872331./1001029789.,                      0., & 
                         &                     0.,    0., 0., & !k0
               & 1201146811./1299019798., -1028468189./846180014.,     0.,        0.,     8478235783./508512852., &
                         & 1311729495./1432422823., -10304129995./1701304382., -48777925059./3047939560., &
                         &  15336726248./1032824649., -45442868181./3398467696.,  3065993473./597172653., & 
                         &                     0.,    0., 0., & !k11
               &                      1.,   185892177./718116043.,     0.,        0.,    -3185094517./667107341., &
                         & -477755414./1098053517.,    -703635378./230739211.,   5731566787./1027545527., &
                         &    5232866602./850066563.,   -4093664535./808688257., 3962137247./1805957418., & 
                         &   65686358./487910083.,    0., 0., & !k12
               &                      1.,   403863854./491063109.,     0.,        0.,    -5068492393./434740067., &
                         &  -411421997./543043805.,     652783627./914296604.,   11173962825./925320556., &
                         & -13158990841./6184727034.,   3936647629./1978049680.,  -160528059./685178525., & 
                         & 248638103./1413531060.,    0., 0., & !k13
               &                      0.,    13451932./455176623.,     0.,        0.,                         0., &
                         &                      0.,    -808719846./976000145.,   1757004468./5645159321., &
                         &     656045339./265891186.,  -3867574721./1518517206.,   465885868./322736535., & 
                         &   53011238./667516719., 2/45., 0.  & !y
               & /)


            call solve_embeed (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, reshape (m, shape=(/14,14/)), maux, 7, 8, ynew)
        end subroutine Dormand_Prince7_8

        ! RK Adap

        recursive subroutine rk_adap_caller (t, y, dt_adap, dydt, integ, p, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            integer, intent(in)                      :: p
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            procedure(integ_tem)                     :: integ
            real*8, dimension(:), intent(in)         :: y
            real*8, dimension(size (y))              :: yaux
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, intent(inout)                    :: dt_adap, dt_used
            real*8                                   :: e_calc

            dt_adap  = max (dt_adap, dt_min)
            call integ (t, y,       dt_adap, dydt, ynew)
            call integ (t, y, 0.5 * dt_adap, dydt, yaux)

            e_calc =  norm2 (ynew - yaux) / (2.**p - 1.)
            if (e_calc < e_tol) then
                dt_used = dt_adap
                dt_adap = max (dt_adap * min (beta * (e_tol / e_calc)**(1./real (p)), 2.), dt_min)
            else
                dt_adap = beta * dt_adap * (e_tol / e_calc)**(1./real (p + 1))
                if ((isnan (dt_adap)) .or. (dt_adap <= dt_min)) then
                    dt_used = dt_min
                    dt_adap = dt_min
                    call integ (t, y, dt_adap, dydt, ynew)
                else
                    call rk_adap_caller (t, y, dt_adap, dydt, integ, p, e_tol, beta, dt_min, dt_used, ynew)
                end if
            end if
        end subroutine rk_adap_caller

        !---------------------------------------------------------------------------------------------
        ! WRAPPERS:
        !---------------------------------------------------------------------------------------------

        !---------------------------------
        !   Call an RK (not embedded nor implicit) integrator
        !---------------------------------

        subroutine integ_caller (t, y, dt, dydt, integ, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            procedure(integ_tem)                     :: integ
            real*8, dimension(size (y)), intent(out) :: ynew

            call integ (t, y, dt, dydt, ynew)
        end subroutine integ_caller

        !---------------------------------
        !   Call an implicit integrator
        !---------------------------------
        
        subroutine implicit_caller (t, y, dt, dydt, integ, max_iter, e_tol, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt, e_tol
            integer, intent(in)                      :: max_iter
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            procedure(implicit_tem)                  :: integ
            real*8, dimension(size (y)), intent(out) :: ynew

            call integ (t, y, dt, dydt, max_iter, e_tol, ynew)
        end subroutine implicit_caller

        !---------------------------------
        !   Call an embedded integrator
        !---------------------------------
        
        subroutine embedded_caller (t, y, dt_adap, dydt, integ, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            real*8, intent(inout)                    :: dt_adap, dt_used
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            procedure(embedded_tem)                  :: integ
            real*8, dimension(size (y)), intent(out) :: ynew

            call integ (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, ynew)
        end subroutine embedded_caller

        !---------------------------------
        !   From =>(f (t, y__), ...)
        !   To: Faux(t, y__)
        !---------------------------------

        ! Here dydt_vec is a pointer to an array of (f__1, ..., f__N);
        !  so this is kind of a wrapper.
        ! dydt_tem_w (t, y__, =>f__)
        !   --> (f__1 (t, y__), ..., f__N (t, y__)) 
        !   --> (der_1, ..., der_N) = der__
        function dydt_wrapper (t, y, dydt) result (der)
            implicit none
            real*8, intent(in)                :: t
            real*8, dimension(:), intent(in)  :: y
            type(dydt_i), dimension(size (y)) :: dydt
            real*8, dimension(size (y))       :: der
            integer                           :: i

            do i = 1, size (y)
                der(i) = dydt(i)%f_i (t , y)
            end do
        end function dydt_wrapper

        ! Here dydt_vec is a pointer to an array of (f__1, ..., f__N)
        !  so this is kind of a wrapper.
        ! integ_tem_w (t, y__, dt, =>f__, integ, ynew__) -->
        !   CREATE Faux (t, y__) = (f__1 (t, y__), ..., f__N (t, y__)) = der__
        !     --> integ (t, y__, dt, Faux__, ynew__) --> ynew__
        subroutine integ_wrapper (t, y, dt, dydt, integ, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            type(dydt_i), dimension(size (y))        :: dydt
            procedure(integ_tem)                     :: integ
            real*8, dimension(size (y)), intent(out) :: ynew

            call integ (t, y, dt, Faux, ynew)

            contains

                function Faux (t, y) result (der)
                    implicit none
                    real*8, intent(in)               :: t
                    real*8, dimension(:), intent(in) :: y
                    real*8, dimension(size (y))      :: der

                    der = dydt_wrapper (t, y, dydt)
                end function Faux
        end subroutine integ_wrapper

        ! Same as before, but for implicit_tem_wrapper integrator
        subroutine implicit_tem_wrapper (t, y, dt, dydt, integ, max_iter, e_tol, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt, e_tol
            integer, intent(in)                      :: max_iter
            real*8, dimension(:), intent(in)         :: y
            type(dydt_i), dimension(size (y))        :: dydt
            procedure(implicit_tem)                  :: integ
            real*8, dimension(size (y)), intent(out) :: ynew

            call integ (t, y, dt, Faux, max_iter, e_tol, ynew)

            contains
            
                function Faux (t, y) result (der)
                    implicit none
                    real*8, intent(in)               :: t
                    real*8, dimension(:), intent(in) :: y
                    real*8, dimension(size (y))      :: der

                    der = dydt_wrapper (t, y, dydt)
                end function Faux
        end subroutine implicit_tem_wrapper

        ! Same as before, but for embedded_tem_wrapper integrator
        subroutine embedded_tem_wrapper (t, y, dt_adap, dydt, integ, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            real*8, intent(inout)                    :: dt_adap, dt_used
            real*8, dimension(:), intent(in)         :: y
            type(dydt_i), dimension(size (y))        :: dydt
            procedure(embedded_tem)                 :: integ
            real*8, dimension(size (y)), intent(out) :: ynew

            call integ (t, y, dt_adap, Faux, e_tol, beta, dt_min, dt_used, ynew)

            contains
            
                function Faux (t, y) result (der)
                    implicit none
                    real*8, intent(in)               :: t
                    real*8, dimension(:), intent(in) :: y
                    real*8, dimension(size (y))      :: der

                    der = dydt_wrapper (t, y, dydt)
                end function Faux
        end subroutine embedded_tem_wrapper

        ! Same as before, but for rk_adap integrator
        subroutine rk_adap_wrapper (t, y, dt_adap, dydt, integ, p, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            integer, intent(in)                      :: p
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            procedure(integ_tem)                     :: integ
            real*8, dimension(:), intent(in)         :: y
            type(dydt_i), dimension(size (y))        :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, intent(inout)                    :: dt_adap, dt_used
            
            call rk_adap_caller (t, y, dt_adap, Faux, integ, p, e_tol, beta, dt_min, dt_used, ynew)

            contains
            
                function Faux (t, y) result (der)
                    implicit none
                    real*8, intent(in)               :: t
                    real*8, dimension(:), intent(in) :: y
                    real*8, dimension(size (y))      :: der

                    der = dydt_wrapper (t, y, dydt)
                end function Faux
        end subroutine rk_adap_wrapper

end module integrators
