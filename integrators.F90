module integrators

    implicit none
    type dydt_i
        procedure(dydt_i_tem), pointer, nopass :: f_i => null ()
    end type dydt_i
    real*8 :: max_iter, e_tol, dt_min, beta, e_calc  ! For adaptive step and implicit

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

        ! in (t, y__, dt, f__, ynew__) -> ynew__, dt
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

        ! Same as before, but for rec_rk4_5
        subroutine rec_rk4_5_w (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, ynew)
            import :: dydt_i
            implicit none
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            real*8, intent(inout)                    :: dt_adap, dt_used 
            real*8, dimension(:), intent(in)         :: y
            type(dydt_i), dimension(size (y))        :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
        end subroutine rec_rk4_5_w
        
        ! Same as before, but for rec_rk_adap
        subroutine rec_rk_adap_w (t, y, dt_adap, dydt, integ, p, e_tol, beta, dt_min, dt_used, ynew)
            import :: dydt_i
            integer, intent(in)                      :: p
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            procedure(integ_tem)                     :: integ
            real*8, dimension(:), intent(in)         :: y
            type(dydt_i), dimension(size (y))        :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, intent(inout)                    :: dt_adap, dt_used        
        end subroutine rec_rk_adap_w

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

        subroutine rksolve_1D (t, y, dt, dydt, m, ynew)
            implicit none
            real*8, intent(in)                   :: t, y, dt
            real*8, dimension(:,:), intent(in)   :: m
            real*8, dimension((size (m, 1) - 1)) :: rk
            procedure(dydt_tem_1D)               :: dydt
            real*8, intent(out)                  :: ynew
            
            call get_rks_1D (t, y, dt, dydt, m, rk)
            ynew = y + dt * dot_product (m(2:,size (m, 1)), rk)
        end subroutine rksolve_1D

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
            
            call rksolve_1D (t, y, dt, dydt, reshape (m, (/7,7/)), ynew)
        end subroutine rungek6_1D

        recursive subroutine rec_rk4_5_1D (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, ynew)
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
!                 dt_adap = max (beta * dt_adap * (e_tol / e_calc)**0.25, dt_min)
            else
                dt_adap = beta * dt_adap * (e_tol / e_calc)**0.2
                if ((isnan (dt_adap)) .or. (dt_adap < dt_min)) then
                    dt_used = dt_min
                    dt_adap = dt_min
                    call rksolve_1D (t, y, dt_adap, dydt, reshape (m, (/7,7/)), ynew)
                else
                    call rec_rk4_5_1D (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, ynew)
                end if 
            end if
        end subroutine rec_rk4_5_1D

        recursive subroutine rec_rk_adap_1D (t, y, dt_adap, dydt, integ, p, e_tol, beta, dt_min, dt_used, ynew)
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
!                 dt_adap = max (beta * dt_adap * (e_tol / e_calc)**(1./real (p)), dt_min)                
            else
                dt_adap = beta * dt_adap * (e_tol / e_calc)**(1./real (p + 1))
                if ((isnan (dt_adap)) .or. (dt_adap <= dt_min)) then
                    dt_used = dt_min
                    dt_adap = dt_min
                    call integ (t, y, dt_adap, dydt, ynew)
                else
                    call rec_rk_adap_1D (t, y, dt_adap, dydt, integ, p, e_tol, beta, dt_min, dt_used, ynew)
                end if
            end if
        end subroutine rec_rk_adap_1D

        !---------------------------------------------------------------------------------------------
        ! ND -> ND
        !---------------------------------------------------------------------------------------------

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
                do j = 1, i           
                    rkaux = rkaux + m(1+j,i) * rk(j,:)                    
                end do
                rk(i, :) = dydt (t + m(1,i) * dt, y + dt * rkaux)
            end do
        end subroutine get_rks

        subroutine rksolve (t, y, dt, dydt, m, ynew)
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
        end subroutine rksolve

        subroutine euler_forward (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(4), parameter          :: m = &
              & (/  0., 0., & !k1
                  & 0., 1.  & !y
              & /)
            
            call rksolve (t, y, dt, dydt, reshape (m, (/2,2/)), ynew)
        end subroutine euler_forward

        subroutine euler_backward (t, y, dt, dydt, max_iter, e_tol, ynew)
            implicit none
            integer, intent(in)                      :: max_iter
            real*8, intent(in)                       :: t, dt, e_tol
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: y1, dy1
            integer                                  :: i
            
            y1 = y
            do i = 1, max_iter
                dy1 = dydt (t + dt, y1) * dt
                y1  = y1  + dy1
                if (norm2 (dy1) <= e_tol) then
                    exit
                end if
            end do
            ynew = y + 0.5 * dydt (t + dt, y1) * dt
        end subroutine euler_backward

        subroutine euler_centred (t, y, dt, dydt, max_iter, e_tol, ynew)
            implicit none
            integer, intent(in)                      :: max_iter
            real*8, intent(in)                       :: t, dt, e_tol
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: y1, dy1
            integer                                  :: i
            
            y1 = y
            do i = 1, max_iter
                dy1 = dydt (t + dt, y1) * dt
                y1  = y1  + dy1
                if (norm2 (dy1) <= e_tol) then
                    exit
                end if
            end do
            ynew = y + 0.5 * (dydt (t + dt, y1) + dydt (t, y)) * dt
        end subroutine euler_centred

        subroutine rungek2 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(9), parameter          :: m = &
              & (/  0.,  0.,  0., & !k1
                  & 1.,  1.,  0., & !k2
                  & 0., 0.5, 0.5  & !y
              & /)
            
            call rksolve (t, y, dt, dydt, reshape (m, (/3,3/)), ynew)
        end subroutine rungek2

        subroutine midpoint2 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(9), parameter          :: m = &
              & (/   0.,  0., 0., & !k1
                  & 0.5, 0.5, 0., & !k2
                  &  0.,  0., 1.  & !y
              & /)
            
            call rksolve (t, y, dt, dydt, reshape (m, (/3,3/)), ynew)
        end subroutine midpoint2

        subroutine strange2 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(9), parameter          :: m = &

              & (/    0.,   0.,   0., & !k1
                  & 0.75, 0.75,   0., & !k2
                  &   0., 1/3., 2/3.  & !y
              & /)
            
            call rksolve (t, y, dt, dydt, reshape (m, (/3,3/)), ynew)
        end subroutine strange2

        subroutine ralston2 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(9), parameter          :: m = &

              & (/    0.,   0.,   0., & !k1
                  & 2/3., 2/3.,   0., & !k2
                  &   0., 0.25, 0.75  & !y
              & /)
            
            call rksolve (t, y, dt, dydt, reshape (m, (/3,3/)), ynew)
        end subroutine ralston2

        subroutine rungek3 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(16), parameter         :: m = &

              & (/    0.,   0.,   0.,   0., & !k1
                  &  0.5,  0.5,   0.,   0., & !k2
                  &   1.,  -1.,   2.,   0., & !k3
                  &   0., 1/6., 2/3., 1/6.  & !y
              & /)
            
            call rksolve (t, y, dt, dydt, reshape (m, (/4,4/)), ynew)
        end subroutine rungek3

        subroutine heun3 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(16), parameter         :: m = &

              & (/    0.,   0.,   0.,   0., & !k1
                  & 1/3., 1/3.,   0.,   0., & !k2
                  & 2/3.,   0., 1/3.,   0., & !k3
                  &   0., 0.25, 1/3., 0.75  & !y
              & /)
            
            call rksolve (t, y, dt, dydt, reshape (m, (/4,4/)), ynew)
        end subroutine heun3

        subroutine ralston3 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(16), parameter         :: m = &

              & (/    0.,   0.,   0.,   0., & !k1
                  &  0.5,  0.5,   0.,   0., & !k2
                  & 0.75,   0., 0.75,   0., & !k3
                  &   0., 2/9., 1/3., 4/9.  & !y
              & /)
            
            call rksolve (t, y, dt, dydt, reshape (m, (/4,4/)), ynew)
        end subroutine ralston3

        subroutine ssprk3 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(16), parameter         :: m = &

              & (/    0.,   0.,   0.,   0., & !k1
                  &   1.,   1.,   0.,   0., & !k2
                  &  0.5, 0.25, 0.25,   0., & !k3
                  &   0., 1/6., 1/6., 2/3.  & !y
              & /)
            
            call rksolve (t, y, dt, dydt, reshape (m, (/4,4/)), ynew)
        end subroutine ssprk3

        subroutine rungek4 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(25), parameter         :: m = &
              & (/    0.,   0.,   0.,   0.,   0., & !k1
                  &  0.5,  0.5,   0.,   0.,   0., & !k2
                  &  0.5,   0.,  0.5,   0.,   0., & !k3
                  &   1.,   0.,   0.,   1.,   0., & !k4
                  &   0., 1/6., 1/3., 1/3., 1/6.  & !y
              & /)
            
            call rksolve (t, y, dt, dydt, reshape (m, (/5,5/)), ynew)
        end subroutine rungek4

        subroutine rungek4_3oct (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(25), parameter         :: m = &

              & (/    0.,    0.,    0.,    0.,    0., & !k1
                  & 1/3.,  1/3.,    0.,    0.,    0., & !k2
                  & 2/3., -1/3.,    1.,    0.,    0., & !k3
                  &   1.,    1.,   -1.,    1.,    0., & !k4
                  &   0., 0.125, 0.375, 0.375, 0.125  & !y
              & /)
            
            call rksolve (t, y, dt, dydt, reshape (m, (/5,5/)), ynew)
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
            
            call rksolve (t, y, dt, dydt, reshape (m, (/7,7/)), ynew)
        end subroutine rungek6

        recursive subroutine rec_rk4_5 (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            real*8, intent(inout)                    :: dt_adap, dt_used
            real*8, dimension(:), intent(in)         :: y
            real*8, dimension(size (y))              :: yaux
            real*8, dimension(6, size (y))           :: rk
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew            
            real*8                                   :: e_calc
            integer                                  :: i
            real*8, parameter, dimension(49)         :: m = &
               & (/  0.,         0.,          0.,          0.,           0.,      0.,    0., & !k1
               &   0.25,       0.25,          0.,          0.,           0.,      0.,    0., & !k2
               &  0.375,      3/32.,       9/32.,          0.,           0.,      0.,    0., & !k3
               & 12/13., 1932/2197., -7200/2197.,  7296/2197.,           0.,      0.,    0., & !k4
               &     1.,   439/216.,         -8.,   3680/513.,   -845/4104.,      0.,    0., & !k5
               &    0.5,     -8/27.,          2., -3544/2565.,   1859/4104., -11/40.,    0., & !k6
               &     0.,    16/135.,          0., 6656/12825., 28561/56430.,   -0.18, 2/55. /) !y

            dt_adap = max (dt_adap, dt_min)
            call get_rks (t, y, dt_adap, dydt, reshape (m, shape=(/7,7/)), rk)
            
            do i = 1, size (y)
                yaux(i) = y(i) + dt_adap * dot_product ((/25/216., 0.,  1408/2565.,   2197/4104.,  -0.2,    0./), rk(:,i))
                ynew(i) = y(i) + dt_adap * dot_product ((/16/135., 0., 6656/12825., 28561/56430., -0.18, 2/55./), rk(:,i))
            end do


            e_calc = norm2 (ynew - yaux)
            if (e_calc < e_tol) then
                dt_used = dt_adap
                dt_adap = max (dt_adap * min (beta * (e_tol / e_calc)**0.25, 2.), dt_min)
                ! dt_adap = max (beta * dt_adap * (e_tol / e_calc)**0.25, dt_min)
                
            else
                dt_adap = beta * dt_adap * (e_tol / e_calc)**0.2
                if ((isnan (dt_adap)) .or. (dt_adap < dt_min)) then
                    dt_adap = dt_min
                    dt_used = dt_min
                    call rksolve (t, y, dt_adap, dydt, reshape (m, (/7,7/)), ynew)
                else
                    call rec_rk4_5 (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, ynew)
                end if 
            end if
        end subroutine rec_rk4_5


        recursive subroutine rec_rk_adap (t, y, dt_adap, dydt, integ, p, e_tol, beta, dt_min, dt_used, ynew)
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
                dt_adap = max (min (beta * dt_adap * (e_tol / e_calc)**(1./real (p)), dt_adap * 2.), dt_min)
                ! dt_adap = max (beta * dt_adap * (e_tol / e_calc)**(1./real (p)), dt_min)
            else
                dt_adap = beta * dt_adap * (e_tol / e_calc)**(1./real (p + 1))
                if ((isnan (dt_adap)) .or. (dt_adap <= dt_min)) then
                    dt_used = dt_min
                    dt_adap = dt_min
                    call integ (t, y, dt_adap, dydt, ynew)
                else
                    call rec_rk_adap (t, y, dt_adap, dydt, integ, p, e_tol, beta, dt_min, dt_used, ynew)
                end if
            end if
        end subroutine rec_rk_adap

        !---------------------------------------------------------------------------------------------
        ! WRAPPERS:
        !---------------------------------------------------------------------------------------------

        !---------------------------------
        !   Call an integrator
        !---------------------------------

        subroutine integ_caller (t, y, dt, dydt, integ, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            procedure(integ_tem)                     :: integ
            real*8, dimension(size (y)), intent(out) :: ynew

            call integ(t, y, dt, dydt, ynew)
        end subroutine integ_caller

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

        ! Same as before, but for rec_rk4_5_wrapper integrator
        subroutine rec_rk4_5_wrapper (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            real*8, intent(inout)                    :: dt_adap, dt_used
            real*8, dimension(:), intent(in)         :: y
            type(dydt_i), dimension(size (y))        :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew

            call rec_rk4_5 (t, y, dt_adap, Faux, e_tol, beta, dt_min, dt_used, ynew)

            contains
            
                function Faux (t, y) result (der)
                    implicit none
                    real*8, intent(in)               :: t
                    real*8, dimension(:), intent(in) :: y
                    real*8, dimension(size (y))      :: der

                    der = dydt_wrapper (t, y, dydt)
                end function Faux

        end subroutine rec_rk4_5_wrapper

        ! Same as before, but for rec_rk_adap integrator
        subroutine rec_rk_adap_wrapper (t, y, dt_adap, dydt, integ, p, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            integer, intent(in)                      :: p
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            procedure(integ_tem)                     :: integ
            real*8, dimension(:), intent(in)         :: y
            type(dydt_i), dimension(size (y))        :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, intent(inout)                    :: dt_adap, dt_used
            
            call rec_rk_adap (t, y, dt_adap, Faux, integ, p, e_tol, beta, dt_min, dt_used, ynew)

            contains
            
                function Faux (t, y) result (der)
                    implicit none
                    real*8, intent(in)               :: t
                    real*8, dimension(:), intent(in) :: y
                    real*8, dimension(size (y))      :: der

                    der = dydt_wrapper (t, y, dydt)
                end function Faux

        end subroutine rec_rk_adap_wrapper

end module integrators