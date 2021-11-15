module integrators

    implicit none
    type dydt_i
        procedure(dydt_i_tem), pointer, nopass :: f_i => null ()
    end type dydt_i
    ! For adaptive step and implicit (might be overwritten)
    integer*4         :: MAX_N_ITER = 5000
    real*8, parameter :: MAX_DT_FAC = 3.
    real*8            :: beta = 0.9, e_tol = 1e-8

    abstract interface

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
        
        ! in (t, y__, dt, f__,  osol, oerr, yaux__, ynew__) -> osol, oerr, yaux__, ynew__
        ! Remember that, in this case,
        !  f__ == (f_1, ..., f_N) must be
        !  pre-defined explicitly
        subroutine embedded_tem (t, y, dt, dydt, osol, oerr, yaux, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            integer*4, intent(out)                   :: osol, oerr
            real*8, dimension(size (y)), intent(out) :: ynew, yaux
        end subroutine embedded_tem

!         ! in (t, y__, dt_adap, f__, e_tol, beta, dt_min, dt_used, ynew) -> ynew__, dt_used
!         ! Remember that, in this case,
!         !  f__ == (f_1, ..., f_N) must be
!         !  pre-defined explicitly
!         subroutine embedded_tem (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, ynew)
!             implicit none
!             real*8, intent(in)                       :: t, e_tol, beta, dt_min
!             real*8, intent(inout)                    :: dt_adap, dt_used
!             real*8, dimension(:), intent(in)         :: y
!             procedure(dydt_tem)                      :: dydt
!             real*8, dimension(size (y)), intent(out) :: ynew
!         end subroutine embedded_tem

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

        ! Same as before, but for embedded_integrators
        subroutine embedded_tem_w (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, integ, ynew)
            import :: dydt_i
            implicit none
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            real*8, intent(inout)                    :: dt_adap, dt_used 
            real*8, dimension(:), intent(in)         :: y
            type(dydt_i), dimension(size (y))        :: dydt
            procedure(embedded_tem)                  :: integ
            real*8, dimension(size (y)), intent(out) :: ynew
        end subroutine embedded_tem_w
        
        ! Same as before, but for rk_adap
        subroutine rk_adap_w (t, y, dt_adap, dydt, integ, p, e_tol, beta, dt_min, dt_used, ynew)
            import :: dydt_i
            integer*4, intent(in)                    :: p
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
        ! ND -> ND
        !---------------------------------------------------------------------------------------------

        !------------------------------------------ SOLVERS ------------------------------------------

        !! Runge Kutta Methods Solver

        subroutine get_rks (t, y, dt, dydt, m, rk)
            implicit none
            real*8, intent(in)                                         :: t, dt
            real*8, dimension(:), intent(in)                           :: y
            real*8, dimension(:,:), intent(in)                         :: m
            real*8, dimension((size (m,1) - 1), size (y)), intent(out) :: rk ! In columns
            real*8, dimension(size (y))                                :: rkaux
            procedure(dydt_tem)                                        :: dydt
            integer*4                                                  :: i, j
            
            do i = 1, size (m, 1) - 1 ! Rows
                rkaux = 0.
                do j = 1, i ! Cols (<i bc its inf triang)    
                    rkaux = rkaux + m(1+j,i) * rk(j,:)                    
                end do
                rk(i, :) = dydt (t + m(1,i) * dt, y + dt * rkaux)
            end do
        end subroutine get_rks
        
        subroutine solve_rk (t, y, dt, dydt, m, ynew)
            implicit none
            real*8, intent(in)                            :: t, dt
            real*8, dimension(:), intent(in)              :: y
            real*8, dimension(:,:), intent(in)            :: m
            real*8, dimension((size (m,1) - 1), size (y)) :: rk
            procedure(dydt_tem)                           :: dydt
            real*8, dimension(size (y)), intent(out)      :: ynew
            integer*4                                     :: i 
            
            call get_rks (t, y, dt, dydt, m, rk)
            
            do i = 1, size (y)
                ynew(i) = y(i) + dt * dot_product ((/m(2:, size (m, 1))/), rk(:,i))
            end do
        end subroutine solve_rk

        !!

        !! Embedded Methods Solver

        !!

        !---------------------------------------- INTEGRATORS ----------------------------------------

        !! Runge Kutta Methods (implicit and explicit)
        
        subroutine Euler1 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            
            ynew = y + dt * dydt (t, y)
        end subroutine Euler1
        
        subroutine Euler_back1 (t, y, dt, dydt, ynew) ! Implicit
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(4), parameter          :: m = &
              & (/1., 1., & !k1
              &   0., 1.  & !y
              & /)
            call solve_rk (t, y, dt, dydt, reshape (m, (/2,2/)), ynew) ! Implicit
        end subroutine Euler_back1

        subroutine Euler2 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(4), parameter          :: m = &
              & (/0.5, 0.5, & !k1
              &    0.,  1.  & !y
              & /)
            
            call solve_rk (t, y, dt, dydt, reshape (m, (/2,2/)), ynew)
        end subroutine Euler2
        
        subroutine Crank_Nicolson2 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: k1, k2
            
            k1 = dydt (t,      y)
            k2 = dydt (t + dt, y + dt * 0.5 * (k1 + k2))
            
            ynew = y + dt * 0.5 * (k1 + k2)
        end subroutine Crank_Nicolson2

        subroutine Heun2 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew              
            real*8, dimension(size (y))              :: k1, k2
            
            k1 = dydt (t,      y)
            k2 = dydt (t + dt, y + dt * k1)
            
            ynew = y + dt * 0.5 * (k1 + k2)
        end subroutine Heun2

        subroutine midpoint2 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: k1, k2
            
            k1 = dydt (t,            y)
            k2 = dydt (t + dt * 0.5, y + dt * 0.5 * k1)
            
            ynew = y + dt * k2
        end subroutine midpoint2

        subroutine strange2 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: k1, k2
            
            k1   = dydt (t,             y)
            k2   = dydt (t + dt * 0.75, y + dt * 0.75 * k1)
            
            ynew = y + dt * (k1 + 2 * k2)/3.
        end subroutine strange2

        subroutine Ralston2 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: k1, k2
            
            k1 = dydt (t,             y)
            k2 = dydt (t + dt * 2/3., y + dt * 2/3. * k1)
            
            ynew = y + dt * 0.25 * (k1 + 3 * k2) 
        end subroutine Ralston2
        
        subroutine Kraaijevanger_Spijker2 (t, y, dt, dydt, ynew) ! Implicit
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(9), parameter          :: m = &

              & (/0.5,  0.5,  0., & !k1
              &   1.5, -0.5,  2., & !k2
              &    0., -0.5, 1.5  & !y
              & /)
            
            call solve_rk (t, y, dt, dydt, reshape (m, (/3,3/)), ynew) 
        end subroutine Kraaijevanger_Spijker2
        
        subroutine Qin_Zhang2 (t, y, dt, dydt, ynew) ! Implicit
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(9), parameter          :: m = &

              & (/0.25, 0.25,   0., & !k1
              &   0.75,  0.5, 0.25, & !k2
              &     0.,  0.5,  0.5  & !y
              & /)
            
            call solve_rk (t, y, dt, dydt, reshape (m, (/3,3/)), ynew)
        end subroutine Qin_Zhang2

        subroutine Runge_Kutta3 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: k1, k2, k3
            
            k1 = dydt (t,            y)
            k2 = dydt (t + dt * 0.5, y + dt * 0.5 * k1)
            k3 = dydt (t + dt,       y + dt * (2. * k2 - k1))
            
            ynew = y + dt * (k1 + 4 * k2 + k3)/6.
        end subroutine Runge_Kutta3

        subroutine Heun3 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: k1, k2, k3
            
            k1 = dydt (t,             y)
            k2 = dydt (t + dt/3.,     y + dt * k1/3.)
            k3 = dydt (t + dt * 2/3., y + dt * k2/3.)
            
            ynew = y + dt * (k1 + 3 * k3) * 0.25
        end subroutine Heun3

        subroutine Ralston3 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: k1, k2, k3

             k1 = dydt (t,             y)
             k2 = dydt (t + dt * 0.5,  y + dt * 0.5 * k1)
             k3 = dydt (t + dt * 0.75, y + dt * 0.75 * k2)
            
             ynew = y + dt * (2 * k1 + 3 * k2 + 4 * k3)/9.
        end subroutine Ralston3

        subroutine SSPRK3 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: k1, k2, k3

             k1 = dydt (t,            y)
             k2 = dydt (t + dt,       y + dt * k1)
             k3 = dydt (t + dt * 0.5, y + dt * 0.25 * (k1 + k2))
            
             ynew = y + dt * (k1 + k2 + 4 * k3)/6.
        end subroutine SSPRK3

        subroutine Ralston4 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: k1, k2, k3, k4
            real*8, dimension(5,5), parameter        :: m = reshape (&

              & (/        0.,         0.,          0.,         0.,         0., & !k1
              &          0.4,        0.4,          0.,         0.,         0., & !k2
              &   0.45573725, 0.29697761,  0.15875964,         0.,         0., & !k3
              &           1., 0.21810040, -3.05096516, 3.83286476,         0., & !k4
              &           0., 0.17476028, -0.55148066, 1.20553560, 0.17118478  & !y
              & /), (/5,5/))
              
             k1 = dydt (t,               y)
             k2 = dydt (t + dt * m(1,2), y + dt * m(2,2) * k1)
             k3 = dydt (t + dt * m(1,3), y + dt * m(2,3) * k1 + m(3,3) * k2)
             k4 = dydt (t + dt,          y + dt * m(2,4) * k1 + m(3,4) * k2 + m(4,4) * k3)
            
             ynew = y + dt * (m(2,5) * k1 + m(3,5) * k2 + m(4,5) * k3 + m(5,5) * k4)
        end subroutine Ralston4
        
        subroutine Gauss_Legendre4 (t, y, dt, dydt, ynew) ! Implicit
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, parameter                        :: sq36 = sqrt (3.)/6.
            real*8, dimension(9), parameter          :: m = &
              & (/0.5 - sq36,        0.25,   0.25 - sq36, & !k1
              &   0.5 + sq36, 0.25 + sq36,          0.25, & !k2
              &               0.,     0.5,           0.5  & !y
              & /)
            
            call solve_rk (t, y, dt, dydt, reshape (m, (/3,3/)), ynew)
        end subroutine Gauss_Legendre4

        subroutine Runge_Kutta4 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: k1, k2, k3, k4
              
             k1 = dydt (t,            y)
             k2 = dydt (t + dt * 0.5, y + dt * 0.5 * k1)
             k3 = dydt (t + dt * 0.5, y + dt * 0.5 * k2)
             k4 = dydt (t + dt,       y + dt * k3)
            
             ynew = y + dt * (k1 + 2 * (k2 + k3) + k4)/6.
        end subroutine Runge_Kutta4

        subroutine Runge_Kutta4_3oct (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: k1, k2, k3, k4
              
             k1 = dydt (t,             y)
             k2 = dydt (t + dt/3.,     y + dt * k1/3.)
             k3 = dydt (t + dt * 2/3., y + dt * (k2 - k1/3.))
             k4 = dydt (t + dt,        y + dt * (k1 + k2 - k3))
            
             ynew = y + dt * (k1 + 3 * (k2 + k3) + k4) * 0.125
        end subroutine Runge_Kutta4_3oct
        
        subroutine Gauss_Legendre6 (t, y, dt, dydt, ynew) ! Implicit
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, parameter                        :: sq15 = sqrt (15.)
            real*8, dimension(16), parameter         :: m = &
              & (/0.5 - sq15 * 0.1,            5/36., 2/9. - sq15/15., 5/36. - sq15/30., & !k1
              &                0.5, 5/36. + sq15/24.,            2/9., 5/36. - sq15/24., & !k2
              &   0.5 + sq15 * 0.1, 5/36. + sq15/30., 2/9. + sq15/15.,            5/36., & !k3
              &                 0.,            5/18.,            4/9.,            5/18.  & !y
              & /)
            
            call solve_rk (t, y, dt, dydt, reshape (m, (/4,4/)), ynew)
        end subroutine Gauss_Legendre6

        subroutine Runge_Kutta6 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: k1, k2, k3, k4, k5, k6
              
             k1 = dydt (t,             y)
             k2 = dydt (t + dt * 0.25, y + dt * 0.25 * k1)
             k3 = dydt (t + dt * 0.25, y + dt * 0.125 * (k1 + k2))
             k4 = dydt (t + dt * 0.5,  y + dt * (k3 - 0.5 * k2))
             k5 = dydt (t + dt * 0.75, y + dt * (0.1875 * k1 + 0.5625 * k4))
             k6 = dydt (t + dt,        y + dt * (-3 * k1 + 2 * k2 + 12 * (k3 - k4) + 8 * k5)/7.)
            
             ynew = y + dt * (7 * (k1 + k6) + 32 * (k3 + k5) + 12 * k4)/90.
        end subroutine Runge_Kutta6

        !!

        !! Embedded

        subroutine Heun_Euler2_1 (t, y, dt, dydt, osol, oaux, yaux, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            integer*4, intent(out)                   :: osol, oaux
            real*8, dimension(size (y)), intent(out) :: ynew, yaux
            real*8, dimension(size (y))              :: k1, k2
            
            osol = 2
            oaux = 1
            
            k1 = dydt (t,      y)
            k2 = dydt (t + dt, y + dt * k1)
            
            ynew = y + dt * 0.5 * (k1 + k2)
            yaux = y + dt * k1
        end subroutine Heun_Euler2_1

        subroutine Fehlberg1_2 (t, y, dt, dydt, osol, oaux, yaux, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            integer*4, intent(out)                   :: osol, oaux
            real*8, dimension(size (y)), intent(out) :: ynew, yaux
            real*8, dimension(size (y))              :: k1, k2, k3
            
            osol = 1
            oaux = 2
            
            k1 = dydt (t,            y)
            k2 = dydt (t + dt * 0.5, y + dt * 0.5 * k1)
            k3 = dydt (t + dt,       y + dt * (k1 + 255 * k2)/256.)
            
            ynew = y + dt * (k1 + 255 * k2)/256.
            yaux = y + dt * (k1 + 500 * k2 + k3)/512.
        end subroutine Fehlberg1_2

        subroutine Bogacki_Shampine3_2  (t, y, dt, dydt, osol, oaux, yaux, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            integer*4, intent(out)                   :: osol, oaux
            real*8, dimension(size (y)), intent(out) :: ynew, yaux
            real*8, dimension(size (y))              :: k1, k2, k3, k4
            
            osol = 3
            oaux = 2
            
            k1 = dydt (t,             y)
            k2 = dydt (t + dt * 0.5,  y + dt * 0.5 * k1)
            k3 = dydt (t + dt * 0.75, y + dt * 0.75 * k2)
            k4 = dydt (t + dt,        y + dt * (2 * k1 + 3 * k2 + 4 * k3)/9.)
            
            ynew = y + dt * (2 * k1 + 3 * k2 + 4 * k3)/9.
            yaux = y + dt * (7/24. * k1 + 0.25 * k2 + k3/3. + 0.125 * k4)
        end subroutine Bogacki_Shampine3_2

        subroutine Zonneveld4_3 (t, y, dt, dydt, osol, oaux, yaux, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            integer*4, intent(out)                   :: osol, oaux
            real*8, dimension(size (y)), intent(out) :: ynew, yaux
            real*8, dimension(size (y))              :: k1, k2, k3, k4, k5
            
            osol = 4
            oaux = 3
            
            k1 = dydt (t,             y)
            k2 = dydt (t + dt * 0.5,  y + dt * 0.5 * k1)
            k3 = dydt (t + dt * 0.5,  y + dt * 0.5 * k2)
            k4 = dydt (t + dt,        y + dt * k3)
            k5 = dydt (t + dt * 0.75, y + dt * (5 * k1 + 7 * k2 + 13 * k3 - k4)/32.)
            
            ynew = y + dt * (k1 + 2 * (k2 + k3) + k4)/6.
            yaux = y + dt * (-0.5 * k1 + (7 * (k2 + k3) - 16 * k5)/3. + 13 * k4/16.)
        end subroutine Zonneveld4_3

        subroutine Merson4_5 (t, y, dt, dydt, osol, oaux, yaux, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            integer*4, intent(out)                   :: osol, oaux
            real*8, dimension(size (y)), intent(out) :: ynew, yaux
            real*8, dimension(size (y))              :: k1, k2, k3, k4, k5
            
            osol = 4
            oaux = 5
            
            k1 = dydt (t,            y)
            k2 = dydt (t + dt/3.,    y + dt * k1/3.)
            k3 = dydt (t + dt/3.,    y + dt * (k1 + k2)/6.)
            k4 = dydt (t + dt * 0.5, y + dt * (0.125 * k1 + 0.375 * k3))
            k5 = dydt (t + dt,       y + dt * (0.5 * k1 - 1.5 * k3 + 2 * k4))
            
            ynew = y + dt * (k1 + 4 * k4 + k5)/6.
            yaux = y + dt * (k1 + 3 * k3 + 4 * k4 + 2 * k5) * 0.1
        end subroutine Merson4_5

        subroutine Fehlberg4_5 (t, y, dt, dydt, osol, oaux, yaux, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            integer*4, intent(out)                   :: osol, oaux
            real*8, dimension(size (y)), intent(out) :: ynew, yaux
            real*8, dimension(size (y))              :: k1, k2, k3, k4, k5, k6
            
            osol = 4
            oaux = 5
            
            k1 = dydt (t,               y)
            k2 = dydt (t + dt * 0.25,   y + dt * 0.25 * k1)
            k3 = dydt (t + dt * 0.375,  y + dt * (3 * k1 + 9 * k2)/32.)
            k4 = dydt (t + dt * 12/13., y + dt * (1932 * k1 - 7200 * k2 + 7296 * k3)/2197.)
            k5 = dydt (t + dt,          y + dt * ((8341 * k1 + 29440 * k3 - 845 * k4)/4104. - 8 * k2))
            k6 = dydt (t + dt * 0.5,    y + dt * (-8 * k1/27. + 2 * k2 - 3544 * k3/2565. + 1859 * k4/4104. - 11 * k5/40.))
            
            ynew = y + dt * ((475 * k1 + 2197 * k4)/4104. + 1408 * k3/2565. - 0.2 * k5)
            yaux = y + dt * ((6688 * k1 + 28561 * k4 + 2052 * k6)/56430. + 6656 * k3/12825. - 0.18 * k5)
        end subroutine Fehlberg4_5

        subroutine Cash_Karp5_4 (t, y, dt, dydt, osol, oaux, yaux, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            integer*4, intent(out)                   :: osol, oaux
            real*8, dimension(size (y)), intent(out) :: ynew, yaux
            real*8, dimension(size (y))              :: k1, k2, k3, k4, k5, k6
            
            osol = 5
            oaux = 4
            
            k1 = dydt (t,              y)
            k2 = dydt (t + dt * 0.2,   y + dt * 0.2 * k1)
            k3 = dydt (t + dt * 0.3,   y + dt * 0.025 * (3 * k1 + 9 * k2))
            k4 = dydt (t + dt * 0.6,   y + dt * 0.1 * (3 * k1 - 9 * k2 + 12 * k3))
            k5 = dydt (t + dt,         y + dt * (-11 * k1 + 135 * k2 - 140 * k3 + 70 * k4)/54.)
            k6 = dydt (t + dt * 0.875, y + dt * (3262 * k1 + 37800 * k2 + 4600 * k3 + 44275 * k4 + 6831 * k5))
            
            ynew = y + dt * (37 * k1/378. + 250 * k3/621. + 125 * k4/594. + 512 * k6/1771.)
            yaux = y + dt * ((5650 * k1 + 13525 * k4)/55296. + 18575 * k3/48384. + 277 * k5/14336. + 0.25 * k6)
        end subroutine Cash_Karp5_4

        subroutine Dormand_Prince5_4 (t, y, dt, dydt, osol, oaux, yaux, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            integer*4, intent(out)                   :: osol, oaux
            real*8, dimension(size (y)), intent(out) :: ynew, yaux
            real*8, dimension(size (y))              :: k1, k2, k3, k4, k5, k6, k7
            
            osol = 5
            oaux = 4

            k1 = dydt (t,             y)
            k2 = dydt (t + dt * 0.2,  y + dt * 0.2 * k1)
            k3 = dydt (t + dt * 0.3,  y + dt * 0.025 * (3 * k1 + 9 * k2))
            k4 = dydt (t + dt * 0.8,  y + dt * (44 * k1 - 168 * k2 + 160 * k3)/45.)
            k5 = dydt (t + dt * 8/9., y + dt * (19372 * k1 - 76080 * k2 + 64448 * k3 - 1908 * k4)/6561.)
            k6 = dydt (t + dt,        y + dt * ((9017 * k1 - 34080 * k2 + 882 * k4)/3168. + 46732 * k3/5247. - 5103 * k5/18656.))
            
            ynew = y + dt * ((35 * k1 + 250 * k4)/384. + 500 * k3/1113. - 2187 * k5/6784. + 11 * k6/84.)
            
            k7 = dydt (t + dt, ynew)
            
            yaux = y + dt * ((5179 * k1 + 35370 * k4)/57600. + 7571 * k3/16695. - 92097 * k5/339200. + 187 * k6/2100. + 0.025 * k7)
        end subroutine Dormand_Prince5_4

        subroutine Verner6_5 (t, y, dt, dydt, osol, oaux, yaux, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            integer*4, intent(out)                   :: osol, oaux
            real*8, dimension(size (y)), intent(out) :: ynew, yaux
            real*8, dimension(size (y))              :: k1, k2, k3, k4, k5, k6, k7, k8
            
            osol = 6
            oaux = 5
            
            k1 = dydt (t,              y)
            k2 = dydt (t + dt/6.,      y + dt * k1/6.)
            k3 = dydt (t + dt * 4/15., y + dt * (4 * k1 + 16 * k2)/75.)
            k4 = dydt (t + dt * 2/3.,  y + dt * ((5 * k1 - 16 * k2)/6. + 2.5 * k3))
            k5 = dydt (t + dt * 5/6.,  y + dt * ((-165 * k1 - 425 * k3)/64. + (880 * k2 + 85 * k4)/96.))
            k6 = dydt (t + dt,         y + dt * ((612 * k1 + 88 * k5)/255. - 8 * k2 + (4015 * k3 - 187 * k4)))
            k7 = dydt (t + dt * 1/15., y + dt * ((-8263 * k1 + 24800 * k2)/15000. - 643 * k3/680. - 0.324 * k4 + 2484 * k5/10625.))
            k8 = dydt (t + dt, y + dt * (3501 * k1/1720. + (297275 * k3 - 367200 * k2)/52632. - 319 * k4/2322. + 24068 * k5/84065. &
                  & + 3850 * k7/26703.))
            
            ynew = y + dt * (0.075 * k1 + 875 * k3/2244. + (3703 * k4 + 125 * k7)/11592. + 264 * k5/1955. + 43 * k8/616.)
            yaux = y + dt * ((13 * k1 + 50 * k4)/160. + (2375 * k3 + 408 * k6)/5984. + 12 * k5/85.)
        end subroutine Verner6_5

!         subroutine Fehlberg7_8 (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, ynew)
!             implicit none
!             real*8, intent(in)                       :: t, e_tol, beta, dt_min
!             real*8, intent(inout)                    :: dt_adap, dt_used
!             real*8, dimension(:), intent(in)         :: y
!             procedure(dydt_tem)                      :: dydt
!             real*8, dimension(size (y)), intent(out) :: ynew
!             real*8, parameter, dimension(13)         :: maux = &
!                & (/0., 0., 0., 0., 0., 34/105., 9/35., 9/35., 9/280., 9/280., 0., 41/840., 41/840./)
!             real*8, parameter, dimension(196)        :: m    = &
!                & (/ 0.,          0.,    0.,      0.,        0.,         0.,       0.,         0.,     0.,      0., &
!                          &     0.,      0., 0., 0., & !k1
!                & 2/27.,       2/27.,    0.,      0.,        0.,         0.,       0.,         0.,     0.,      0., &
!                          &     0.,      0., 0., 0., & !k2
!                &  1/9.,       1/36., 1/12.,      0.,        0.,         0.,       0.,         0.,     0.,      0., &
!                          &     0.,      0., 0., 0., & !k3
!                &  1/6.,       1/24.,    0.,   0.125,        0.,         0.,       0.,         0.,     0.,      0., &
!                          &     0.,      0., 0., 0., & !k4
!                & 5/12.,       5/12.,    0., -1.5625,    1.5625,         0.,       0.,         0.,     0.,      0., &
!                          &     0.,      0., 0., 0., & !k5
!                &   0.5,        0.05,    0.,      0.,      0.25,        0.2,       0.,         0.,     0.,      0., &
!                          &     0.,      0., 0., 0., & !k6
!                &  5/6.,    -25/108.,    0.,      0.,  125/108.,    -65/27.,  125/54.,         0.,     0.,      0., &
!                          &     0.,      0., 0., 0., & !k7
!                &  1/6.,     31/300.,    0.,      0.,        0.,    61/225.,    -2/9.,    13/900.,     0.,      0., &
!                          &     0.,      0., 0., 0., & !k8
!                &  2/3.,          2.,    0.,      0.,    -53/6.,    704/45.,  -107/9.,     67/90.,     3.,      0., &
!                          &     0.,      0., 0., 0., & !k9
!                &  1/3.,    -91/108.,    0.,      0.,   23/108.,  -976/135.,  311/54.,    -19/60.,  17/6.,  -1/12., &
!                          &     0.,      0., 0., 0., & !k10
!                &    1.,  2383/4100.,    0.,      0., -341/164., 4496/1025., -301/82., 2133/4100., 45/82., 45/164., &
!                          & 18/41.,      0., 0., 0., & !k11
!                &    0.,      3/205.,    0.,      0.,        0.,         0.,   -6/41.,    -3/205., -3/41.,   3/41., &
!                          &  6/41.,      0., 0., 0., & !k12
!                &    1., -1777/4100.,    0.,      0., -341/164., 4496/1025., -289/82., 2193/4100., 51/82., 33/164., &
!                          & 19/41.,      0., 1., 0., & !k13
!                &    0.,     41/840.,    0.,      0.,        0.,         0.,  34/105.,      9/35.,  9/35.,  9/280., &
!                          & 9/280., 41/840., 0., 0. & !y
!                & /)
! 
!             call solve_embed (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, reshape (m, shape=(/14,14/)), maux, 7, 8, ynew)
!         end subroutine Fehlberg7_8
!         
!         subroutine Dormand_Prince8_7 (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, ynew)
!             implicit none
!             real*8, intent(in)                       :: t, e_tol, beta, dt_min
!             real*8, intent(inout)                    :: dt_adap, dt_used
!             real*8, dimension(:), intent(in)         :: y
!             procedure(dydt_tem)                      :: dydt
!             real*8, dimension(size (y)), intent(out) :: ynew
!             real*8, parameter, dimension(13)         :: maux = &
!                & (/14005451./335480064., 0., 0., 0., 0., -59238493./1068277825., 181606767./758867731., 561292985./797845732.,&
!                          & -1041891430./1371343529., 760417239./1151165299., 118820643./751138087., -528747749./2220607170., 0.25/)
!             real*8, parameter, dimension(196)        :: m    = &
!                & (/                   0.,                      0.,     0.,        0.,                         0., &
!                          &                      0.,                        0.,                        0., &
!                          &                        0.,                        0.,                      0., & 
!                          &                     0.,    0., 0., & !k1
!                &                   1/18.,                   1/18.,     0.,        0.,                         0., &
!                          &                      0.,                        0.,                        0., &
!                          &                        0.,                        0.,                      0., & 
!                          &                     0.,    0., 0., & !k2
!                &                   1/12.,                   1/48., 0.0625,        0.,                         0., &
!                          &                      0.,                        0.,                        0., &
!                          &                        0.,                        0.,                      0., & 
!                          &                     0.,    0., 0., & !k3
!                &                   0.125,                 0.03125,     0.,   0.09375,                         0., &
!                          &                      0.,                        0.,                        0., &
!                          &                        0.,                        0.,                      0., & 
!                          &                     0.,    0., 0., & !k4
!                &                  0.3125,                  0.3125,     0., -1.171875,                   1.171875, &
!                          &                      0.,                        0.,                        0., &
!                          &                        0.,                        0.,                      0., & 
!                          &                     0.,    0., 0., & !k5
!                &                   0.375,                  0.0375,     0.,        0.,                     0.1875, &
!                          &                    0.15,                        0.,                        0., &
!                          &                        0.,                        0.,                      0., & 
!                          &                     0.,    0., 0., & !k6
!                &                  0.1475,    29443841./614563906.,     0.,        0.,       77736538./692538347., &
!                          &  -28693883./1125000000.,     23124283./1800000000.,                        0., &
!                          &                        0.,                        0.,                      0., & 
!                          &                     0.,    0., 0., & !k7
!                &                   0.465,    16016141./946692911.,     0.,        0.,       61564180./158732637., &
!                          &    22789713./633445777.,    545815736./2771057229.,   -180193667./1043307555., &
!                          &                        0.,                        0.,                      0., & 
!                          &                     0.,    0., 0., & !k8
!                & 5490023248./9719169821.,    39632708./573591083.,     0.,        0.,     -433636366./683701615., &
!                          & -421739975./2616292301.,     100302831./723423059.,     790204164./839813087., &
!                          &    800635310./3783071287.,                        0.,                      0., & 
!                          &                     0.,    0., 0., & !k9
!                &                    0.65,  246121993./1340847787.,     0.,        0., -37695042795./15268766246., &
!                          & -309121744./1061227803.,     -12992083./490766935.,   6005943493./2108947869., &
!                          &    393006217./1396673457.,    123872331./1001029789.,                      0., & 
!                          &                     0.,    0., 0., & !k0
!                & 1201146811./1299019798., -1028468189./846180014.,     0.,        0.,     8478235783./508512852., &
!                          & 1311729495./1432422823., -10304129995./1701304382., -48777925059./3047939560., &
!                          &  15336726248./1032824649., -45442868181./3398467696.,  3065993473./597172653., & 
!                          &                     0.,    0., 0., & !k11
!                &                      1.,   185892177./718116043.,     0.,        0.,    -3185094517./667107341., &
!                          & -477755414./1098053517.,    -703635378./230739211.,   5731566787./1027545527., &
!                          &    5232866602./850066563.,   -4093664535./808688257., 3962137247./1805957418., & 
!                          &   65686358./487910083.,    0., 0., & !k12
!                &                      1.,   403863854./491063109.,     0.,        0.,    -5068492393./434740067., &
!                          &  -411421997./543043805.,     652783627./914296604.,   11173962825./925320556., &
!                          & -13158990841./6184727034.,   3936647629./1978049680.,  -160528059./685178525., & 
!                          & 248638103./1413531060.,    0., 0., & !k13
!                &                      0.,    13451932./455176623.,     0.,        0.,                         0., &
!                          &                      0.,    -808719846./976000145.,   1757004468./5645159321., &
!                          &     656045339./265891186.,  -3867574721./1518517206.,   465885868./322736535., & 
!                          &   53011238./667516719., 2/45., 0.  & !y
!                & /)
! 
! 
!             call solve_embed (t, y, dt_adap, dydt, e_tol, beta, dt_min, dt_used, reshape (m, shape=(/14,14/)), maux, 8, 7, ynew)
!         end subroutine Dormand_Prince8_7

        !!

        !! Runge Kutta Half_Step

        recursive subroutine rk_half_step (t, y, dt_adap, dydt, integ, p, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            integer*4, intent(in)                    :: p
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            procedure(integ_tem)                     :: integ
            real*8, dimension(:), intent(in)         :: y
            real*8, dimension(size (y))              :: yhalf, yaux
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, intent(inout)                    :: dt_adap, dt_used
            real*8                                   :: e_calc, ratio, hdt_adap
            integer*4                                :: iter = 1

            dt_adap  = max (dt_adap, dt_min)
            hdt_adap = 0.5 * dt_adap

            call integ (           t,     y,  dt_adap, dydt,  ynew)
            call integ (           t,     y, hdt_adap, dydt, yhalf)                       
            call integ (t + hdt_adap, yhalf, hdt_adap, dydt,  yaux)

            e_calc =  norm2 (ynew - yaux) / real(2**p - 1, kind=8)
            ratio  = e_tol / e_calc
            if ((ratio > 1.) .or. (iter .eq. MAX_N_ITER)) then
                dt_used = dt_adap
                dt_adap = dt_adap * min (beta * ratio**(1. / real (p + 1, kind=8)), MAX_DT_FAC)
                iter    = 1
            else
                dt_adap = dt_adap * min (beta * ratio**(1. / real (p, kind=8)), MAX_DT_FAC)
                if ((isnan (dt_adap)) .or. (dt_adap .le. dt_min)) then
                    dt_used = dt_min
                    dt_adap = dt_min
                    call integ (t, y, dt_adap, dydt, ynew)
                    iter = 1
                else
                    call rk_half_step (t, y, dt_adap, dydt, integ, p, e_tol, beta, dt_min, dt_used, ynew)
                    iter = iter + 1
                end if
            end if
        end subroutine rk_half_step

        !!

        !---------------------------------------------------------------------------------------------
        ! CALLERS:
        !---------------------------------------------------------------------------------------------

        !---------------------------------
        !   Call an RK (not embedded nor implicit) integrator
        !---------------------------------

        subroutine integ_caller (t, y, dt_min, dydt, integ, dt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt, dt_min
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            procedure(integ_tem)                     :: integ
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: yaux
            real*8                                   :: time, t_end, dt_used

            ynew    = y
            time    = t
            t_end   = time + dt
            dt_used = min (dt_min, dt)
            do while (time < t_end)
                yaux = ynew
                call integ (time, yaux, dt_used, dydt, ynew)
                time = time + dt_used
            end do
        end subroutine integ_caller

        !---------------------------------
        !   Call an embedded integrator
        !---------------------------------
        
         recursive subroutine solve_embed (t, y, dt_adap, dydt, integ, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            real*8, intent(in)                            :: t, e_tol, beta, dt_min
            real*8, intent(inout)                         :: dt_adap, dt_used
            real*8, dimension(:), intent(in)              :: y
            integer*4, save                               :: osol, oaux, iter = 1
            real*8, dimension(size (y))                   :: yaux
            procedure(dydt_tem)                           :: dydt
            procedure(embedded_tem)                       :: integ
            real*8, dimension(size (y)), intent(out)      :: ynew            
            real*8                                        :: e_calc, ratio

            
            dt_adap = max (dt_adap, dt_min)
            
            call integ (t, y, dt_adap, dydt, osol, oaux, yaux, ynew)

            e_calc = norm2 (ynew - yaux)
            ratio = e_tol / e_calc
            if ((ratio > 1.) .or. (iter .eq. MAX_N_ITER)) then
                dt_used = dt_adap
                dt_adap = dt_adap * min (beta * ratio**(1. / osol), MAX_DT_FAC)
                iter = 1
            else
                if (dt_adap .eq. dt_min) then
                    dt_used = dt_min
                else
                    dt_adap = dt_adap * min (beta * ratio**(1. / oaux), MAX_DT_FAC)
                    if ((isnan (dt_adap)) .or. (dt_adap < dt_min)) then
                        dt_adap = dt_min
                        dt_used = dt_min
                        call integ (t, y, dt_adap, dydt, osol, oaux, yaux, ynew)
                        iter = 1
                    else
                        call solve_embed (t, y, dt_adap, dydt, integ, e_tol, beta, dt_min, dt_used, ynew)
                        iter = iter + 1
                    end if 
                end if
            end if
        end subroutine solve_embed
        
        subroutine embedded_caller (t, y, dt_adap, dydt, integ, e_tol, beta, dt_min, dt, ynew)
            implicit none
            real*8, intent(in)                       :: t, e_tol, beta, dt_min, dt
            real*8, dimension(:), intent(in)         :: y
            real*8, intent(inout)                    :: dt_adap
            procedure(dydt_tem)                      :: dydt
            procedure(embedded_tem)                  :: integ
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: yaux
            real*8                                   :: time, t_end, dtmin, dtused

            ynew  = y
            time  = t
            t_end = time + dt
            dtmin = min (dt_min, dt)
            do while (time < t_end)
                yaux  = ynew
                dt_adap = min (dt_adap, t_end - time)
                call solve_embed (time, yaux, dt_adap, dydt, integ, e_tol, beta, dtmin, dtused, ynew)
                time = time + dtused                
            end do
        end subroutine embedded_caller

        !---------------------------------
        !   Call rk_half_step integrator
        !---------------------------------

        subroutine rk_half_step_caller (t, y, dt_adap, dydt, integ, p, e_tol, beta, dt_min, dt, ynew)
            implicit none
            real*8, intent(in)                       :: t, e_tol, beta, dt, dt_min
            real*8, dimension(:), intent(in)         :: y
            real*8, intent(inout)                    :: dt_adap
            procedure(dydt_tem)                      :: dydt
            procedure(integ_tem)                     :: integ
            integer*4, intent(in)                    :: p
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: yaux
            real*8                                   :: time, t_end, dtmin, dtused

            ynew  = y
            time  = t
            t_end = time + dt
            dtmin = min (dt_min, dt)
            do while (time < t_end)
                yaux  = ynew
                dt_adap = min (dt_adap, t_end - time)
                call rk_half_step (time, yaux, dt_adap, dydt, integ, p, e_tol, beta, dtmin, dtused, ynew)
                time = time + dtused
            end do
        end subroutine rk_half_step_caller

        !---------------------------------------------------------------------------------------------
        ! WRAPPERS:
        !---------------------------------------------------------------------------------------------

        !---------------------------------
        !   From =>(f (t, y__), ...)
        !   To: Faux(t, y__)
        !---------------------------------

        ! Here dydt is a pointer to an array of (f__1, ..., f__N);
        !  so this is kind of a wrapper.
        ! dydt_w (t, y__, =>f__)
        !   --> (f__1 (t, y__), ..., f__N (t, y__)) 
        !   --> (der_1, ..., der_N) = der__
        function dydt_wrapper (t, y, dydt) result (der)
            implicit none
            real*8, intent(in)                :: t
            real*8, dimension(:), intent(in)  :: y
            type(dydt_i), dimension(size (y)) :: dydt
            real*8, dimension(size (y))       :: der
            integer*4                         :: i

            do i = 1, size (y)
                der(i) = dydt(i)%f_i (t , y)
            end do
        end function dydt_wrapper

        ! Here dydt is a pointer to an array of (f__1, ..., f__N)
        !  so this is kind of a wrapper.
        ! integ_wrapper (t, y__, dt, =>f__, integ, dt_min, ynew__) -->
        !   CREATE Faux (t, y__) = (f__1 (t, y__), ..., f__N (t, y__)) = der__
        !     --> implicit_caller (t, y__, dt, Faux__, integ, dt_min, ynew__) --> ynew__
        subroutine integ_wrapper (t, y, dt, dydt, integ, dt_min, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt, dt_min
            real*8, dimension(:), intent(in)         :: y
            type(dydt_i), dimension(size (y))        :: dydt
            procedure(integ_tem)                     :: integ
            real*8, dimension(size (y)), intent(out) :: ynew

            call integ_caller (t, y, dt, Faux, integ, dt_min, ynew)

            contains

                function Faux (ti, yi) result (der)
                    implicit none
                    real*8, intent(in)               :: ti
                    real*8, dimension(:), intent(in) :: yi
                    real*8, dimension(size (yi))     :: der

                    der = dydt_wrapper (ti, yi, dydt)
                end function Faux
        end subroutine integ_wrapper
        
        ! Same as before, but for embedded_caller
        subroutine embedded_tem_wrapper (t, y, dt_adap, dydt, integ, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            real*8, intent(inout)                    :: dt_adap, dt_used
            real*8, dimension(:), intent(in)         :: y
            type(dydt_i), dimension(size (y))        :: dydt
            procedure(embedded_tem)                  :: integ
            real*8, dimension(size (y)), intent(out) :: ynew

            call embedded_caller (t, y, dt_adap, Faux, integ, e_tol, beta, dt_min, dt_used, ynew)

            contains
            
                function Faux (ti, yi) result (der)
                    implicit none
                    real*8, intent(in)               :: ti
                    real*8, dimension(:), intent(in) :: yi
                    real*8, dimension(size (yi))     :: der

                    der = dydt_wrapper (ti, yi, dydt)
                end function Faux
        end subroutine embedded_tem_wrapper

        ! Same as before, but for rk_half_step_caller
        subroutine rk_adap_wrapper (t, y, dt_adap, dydt, integ, p, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            integer*4, intent(in)                    :: p
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            procedure(integ_tem)                     :: integ
            real*8, dimension(:), intent(in)         :: y
            type(dydt_i), dimension(size (y))        :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, intent(inout)                    :: dt_adap, dt_used
            
            call rk_half_step_caller (t, y, dt_adap, Faux, integ, p, e_tol, beta, dt_min, dt_used, ynew)

            contains
            
                function Faux (ti, yi) result (der)
                    implicit none
                    real*8, intent(in)               :: ti
                    real*8, dimension(:), intent(in) :: yi
                    real*8, dimension(size (yi))     :: der

                    der = dydt_wrapper (ti, yi, dydt)
                end function Faux
        end subroutine rk_adap_wrapper

end module integrators
