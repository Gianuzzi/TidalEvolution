module integrators

    implicit none
    type dydt_i
        procedure(dydt_i_tem), pointer, nopass :: f_i => null ()
    end type dydt_i
    ! For adaptive step and implicit (might be overwritten)
    integer*4         :: MAX_N_ITER = 100
    real*8, parameter :: MAX_DT_FAC = 5., SAFE_LOW = 1e-300
    real*8            :: BETA = 0.9, E_TOL = 1e-8

    ! Aux Constants
    real*8, parameter :: C1_3   = 1/3., C2_3 = C1_3 * 2
    real*8, parameter :: C1_5   = 1/5.
    real*8, parameter :: C1_6   = 1/6., C5_6 = C1_6 * 5
    real*8, parameter :: C1_7   = 1/7.
    real*8, parameter :: C1_9   = 1/9., C2_9 = C1_9 * 2, C4_9 = C1_9 * 4, C8_9 = C1_9 * 8
    real*8, parameter :: C1_12  = 1/12., C5_12 = C1_12 * 5
    real*8, parameter :: C1_15  = 1/15., C4_15 = C1_15 * 4
    real*8, parameter :: C1_18  = 1/18., C5_18 = C1_18 * 5
    real*8, parameter :: C2_27  = 2/27.
    real*8, parameter :: C5_36  = 5/36.
    real*8, parameter :: C2_45  = 2/45.
    real*8, parameter :: C1_48  = 1/48.
    real*8, parameter :: C1_840 = 1/840., C41_840 = C1_840 * 41
    real*8, parameter :: SQ3    = sqrt(3.)
    real*8, parameter :: SQ15   = sqrt(15.)
    real*8, parameter :: SQ3_6  = SQ3 * C1_6, SQ3_2 = SQ3_6 * 2
    real*8, parameter :: SQ15_5 = SQ15 * 0.2, SQ15_10 = SQ15_5 * 2, SQ15_24 = SQ15/24.
    real*8, parameter :: SQ15_30 = SQ15/30., SQ15_15 = SQ15_30 * 2

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
        subroutine embedded_tem (t, y, dt, dydt, osol, oaux, yaux, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            integer*4, intent(out)                   :: osol, oaux
            real*8, dimension(size (y)), intent(out) :: ynew, yaux
        end subroutine embedded_tem

        subroutine calc_rk (t, y, dt, dydt, kin, kout)
            implicit none
            real*8, intent(in)                                          :: t, dt
            procedure(dydt_tem)                                         :: dydt
            real*8, dimension(:, :), intent(in)                         :: kin
            real*8, dimension(size (kin, 1)), intent(in)                :: y
            real*8, dimension(size (kin, 1), size (kin,2)), intent(out) :: kout
        end subroutine calc_rk 

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

        !! Implicit Methods Solver

        subroutine solve_1k_implicit (t, y, dt, dydt, kprev, cte, k)
            implicit none
            real*8, intent(in)                       :: t, dt, cte
            real*8, dimension(:), intent(in)         :: y
            real*8, dimension(size (y)), intent(in)  :: kprev
            real*8, dimension(size (y)), intent(out) :: k
            real*8, dimension(size (y))              :: kaux
            procedure(dydt_tem)                      :: dydt
            integer*4                                :: i
            
            k = dydt (t, y + dt * kprev)
            do i = 1, MAX_N_ITER
                kaux = k
                k    = dydt (t, y + dt * (kprev + cte * k))
                if (maxval( abs ((kaux - k) / (kaux + SAFE_LOW))) .le. E_TOL) then                    
                    exit
                end if
            end do
        end subroutine solve_1k_implicit

        subroutine solve_rk_implicit (t, y, dt, dydt, solver, rk)
            implicit none
            real*8, intent(in)                          :: t, dt
            procedure(dydt_tem)                         :: dydt
            real*8, dimension(:, :), intent(inout)      :: rk
            real*8, dimension(size (rk,1), size (rk,2)) :: rkold
            real*8, dimension(size (rk,1)), intent(in)  :: y
            procedure(calc_rk)                          :: solver
            integer*4                                   :: i

            do i = 1, MAX_N_ITER
                rkold = rk
                call solver (t, y, dt, dydt, rkold, rk)
                if (maxval( abs ((rkold - rk) / (rkold + SAFE_LOW))) .le. E_TOL) then
                    exit
                end if
            end do
        end subroutine solve_rk_implicit

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

        subroutine solve_rk_embed (t, y, dt, dydt, m, maux, yaux, ynew)
            implicit none
            real*8, intent(in)                              :: t, dt
            real*8, dimension(:), intent(in)                :: y
            real*8, dimension(:,:), intent(in)              :: m
            real*8, dimension((size (m,1) - 1), size (y))   :: rk
            real*8, dimension((size (m,1) - 1)), intent(in) :: maux
            procedure(dydt_tem)                             :: dydt
            real*8, dimension(size (y)), intent(out)        :: yaux, ynew
            integer*4                                       :: i 
            
            call get_rks (t, y, dt, dydt, m, rk)
            
            do i = 1, size (y)
                yaux(i) = y(i) + dt * dot_product (maux, rk(:,i))
                ynew(i) = y(i) + dt * dot_product ((/m(2:, size (m, 1))/), rk(:,i))
            end do
        end subroutine solve_rk_embed

        !!

        !! Embedded Methods Solver

        recursive subroutine solve_embed (t, y, dt_adap, dydt, integ, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            real*8, intent(in)                            :: t, e_tol, beta, dt_min
            real*8, intent(inout)                         :: dt_adap, dt_used
            real*8, dimension(:), intent(in)              :: y
            integer*4, save                               :: osol, oaux, iter = 0
            real*8, dimension(size (y))                   :: yaux, yscal
            procedure(dydt_tem)                           :: dydt
            procedure(embedded_tem)                       :: integ
            real*8, dimension(size (y)), intent(out)      :: ynew            
            real*8                                        :: e_calc, ratio

            iter    = iter + 1
            dt_adap = max (dt_adap, dt_min)            
            
            call integ (t, y, dt_adap, dydt, osol, oaux, yaux, ynew)

            yscal = abs (y + dt_adap * dydt (t, y)) + SAFE_LOW
            
            e_calc = max (maxval (abs ((ynew - yaux) / yscal )), SAFE_LOW)
            ratio  = e_tol / e_calc
            if (ratio > 1.) then
                dt_used = dt_adap
                dt_adap = dt_adap * min (beta * ratio**(1. / osol), MAX_DT_FAC)
                iter = 0
            else
                if (dt_adap .eq. dt_min) then
                    dt_used = dt_min
                    iter = 0
                else
                    dt_adap = dt_adap * min (beta * ratio**(1. / oaux), MAX_DT_FAC)
                    if ((isnan (dt_adap)) .or. (dt_adap < dt_min) .or. (iter .eq. MAX_N_ITER)) then
                        dt_adap = dt_min
                        dt_used = dt_min
                        call integ (t, y, dt_adap, dydt, osol, oaux, yaux, ynew)
                        iter = 0
                    else
                        call solve_embed (t, y, dt_adap, dydt, integ, e_tol, beta, dt_min, dt_used, ynew)
                    end if 
                end if
            end if
        end subroutine solve_embed

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
            real*8, dimension(size (y))              :: kaux, k1
            real*8                                   :: aux

            kaux = 0.
            aux  = 1.
            call solve_1k_implicit (t + dt, y, dt, dydt, kaux, aux, k1)

            ynew = y + dt * k1
        end subroutine Euler_back1

        subroutine Euler_center2 (t, y, dt, dydt, ynew) ! Implicit
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: kaux, k1
            real*8                                   :: aux

            kaux = 0.
            aux  = 0.5
            call solve_1k_implicit (t + dt * 0.5, y, dt, dydt, kaux, aux, k1)

            ynew = y + dt * k1
        end subroutine Euler_center2
        
        subroutine Crank_Nicolson2 (t, y, dt, dydt, ynew) ! Implicit
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: kaux, k1, k2
            real*8                                   :: aux

            k1   = dydt (t, y)
            kaux = k1 * 0.5
            aux  = 0.5
            call solve_1k_implicit (t + dt, y, dt, dydt, kaux, aux, k2)

            ynew = y + dt * (k1 + k2) * 0.5
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
            
            ynew = y + dt * (k1 + k2) * 0.5
        end subroutine Heun2

        subroutine midpoint2 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: k1, k2
            
            k1 = dydt (t,            y)
            k2 = dydt (t + dt * 0.5, y + dt * k1 * 0.5)
            
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
            k2   = dydt (t + dt * 0.75, y + dt * k1 * 0.75)
            
            ynew = y + dt * (k1 + k2 * 2) * C1_3
        end subroutine strange2

        subroutine Ralston2 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: k1, k2
            
            k1 = dydt (t,             y)
            k2 = dydt (t + dt * C2_3, y + dt * k1 * C2_3)
            
            ynew = y + dt * (k1 + k2 * 3) * 0.25
        end subroutine Ralston2

        subroutine Hammer_Hollingsworth2 (t, y, dt, dydt, ynew) ! Implicit
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: kaux, k1, k2
            real*8                                   :: aux

            k1 = dydt (t, y)
            kaux = k1 * C1_3
            aux  = C1_3
            call solve_1k_implicit (t + dt * C2_3, y, dt, dydt, kaux, aux, k2)

            ynew = y + dt * (k1 + k2 * 3) * 0.25
        end subroutine Hammer_Hollingsworth2
        
        subroutine Kraaijevanger_Spijker2 (t, y, dt, dydt, ynew) ! Implicit
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: kaux, k1, k2
            real*8                                   :: aux

            kaux = 0.
            aux  = 0.5
            call solve_1k_implicit (t + dt * 0.5, y, dt, dydt, kaux, aux, k1)
            kaux = - k1 * 0.5
            aux  = 2.
            call solve_1k_implicit (t + dt * 1.5, y, dt, dydt, kaux, aux, k2)

            ynew = y + dt * (- k1 + k2 * 3) * 0.5
        end subroutine Kraaijevanger_Spijker2
        
        subroutine Qin_Zhang2 (t, y, dt, dydt, ynew) ! Implicit
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: kaux, k1, k2
            real*8, parameter                        :: aux = 0.25

            kaux = 0.
            call solve_1k_implicit (t + dt * 0.25, y, dt, dydt, kaux, aux, k1)
            kaux = k1 * 0.5
            call solve_1k_implicit (t + dt * 0.75, y, dt, dydt, kaux, aux, k2)

            ynew = y + dt * (k1 + k2) * 0.5
        end subroutine Qin_Zhang2

        subroutine Runge_Kutta3 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: k1, k2, k3
            
            k1 = dydt (t,            y)
            k2 = dydt (t + dt * 0.5, y + dt * k1 * 0.5)
            k3 = dydt (t + dt,       y + dt * (- k1 + k2 * 2))
            
            ynew = y + dt * (k1 + k2 * 4 + k3) * C1_6
        end subroutine Runge_Kutta3

        subroutine Heun3 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: k1, k2, k3
            
            k1 = dydt (t,             y)
            k2 = dydt (t + dt * C1_3, y + dt * k1 * C1_3)
            k3 = dydt (t + dt * C2_3, y + dt * k2 * C2_3)
            
            ynew = y + dt * (k1 + k3 * 3) * 0.25
        end subroutine Heun3

        subroutine Ralston3 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: k1, k2, k3

            k1 = dydt (t,             y)
            k2 = dydt (t + dt * 0.5,  y + dt * k1 * 0.5)
            k3 = dydt (t + dt * 0.75, y + dt * k2 * 0.75)

            ynew = y + dt * (k1 * 2 + k2 * 3 + k3 * 4) * C1_9
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
            k3 = dydt (t + dt * 0.5, y + dt * (k1 + k2) * 0.25)

            ynew = y + dt * (k1 + k2 + k3 * 4) * C1_6
        end subroutine SSPRK3

        subroutine Crouzeix3 (t, y, dt, dydt, ynew) ! Implicit
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: kaux, k1, k2
            real*8                                   :: aux

            kaux = 0.
            aux  = 0.5 + SQ3_6
            call solve_1k_implicit (t + dt * aux, y, dt, dydt, kaux, aux, k1)
            kaux = - k1 * SQ3_6 * 2
            call solve_1k_implicit (t + dt * (0.5 - SQ3_6), y, dt, dydt, kaux, aux, k2)

            ynew = y + dt * (k1 + k2) * 0.5
        end subroutine Crouzeix3

        subroutine Runge_Kutta_implicit3 (t, y, dt, dydt, ynew) ! Implicit
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: kaux, k1, k2, k3, k4
            real*8, parameter                        :: aux = 0.5

            kaux = 0.
            call solve_1k_implicit (t + dt * 0.5, y, dt, dydt, kaux, aux, k1)
            kaux = k1 * C1_6
            call solve_1k_implicit (t + dt * C2_3, y, dt, dydt, kaux, aux, k2)
            kaux = (- k1 + k2) * 0.5
            call solve_1k_implicit (t + dt * 0.5, y, dt, dydt, kaux, aux, k3)
            kaux = ((k1 - k2) * 3 + k3) * 0.5
            call solve_1k_implicit (t + dt, y, dt, dydt, kaux, aux, k4)

            ynew = y + dt * (kaux + k4 * 0.5)
        end subroutine Runge_Kutta_implicit3

        subroutine Ralston4 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: k1, k2, k3, k4

            k1 = dydt (t,                   y)
            k2 = dydt (t + dt * 0.4,        y + dt * k1 * 0.4)
            k3 = dydt (t + dt * 0.45573725, y + dt * (k1 * 0.29697761 + k2 * 0.15875964))
            k4 = dydt (t + dt,              y + dt * (k1 * 0.21810040 - k2 * 3.05096516 + k3 * 3.83286476))

            ynew = y + dt * (k1 * 0.17476028 - k2 * 0.55148066 + k3 * 1.20553560 + k4 * 0.17118478)
        end subroutine Ralston4

        subroutine Lobatto4 (t, y, dt, dydt, ynew) ! Implicit
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: kaux, k1, k2, k3
            real*8, parameter                        :: aux = 0.25

            k1   = dydt (t, y)
            kaux = k1 * 0.25
            call solve_1k_implicit (t + dt * 0.5, y, dt, dydt, kaux, aux, k2)
            k3 = dydt (t + dt, y + dt * k2)

            ynew = y + dt * (k1 + k2 * 4 + k3) * C1_6
        end subroutine Lobatto4

        subroutine Runge_Kutta4 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: k1, k2, k3, k4
            real*8                                   :: aut
              
            aut = dt * 0.5

            k1 = dydt (t,       y)
            k2 = dydt (t + aut, y + dt * k1 * 0.5)
            k3 = dydt (t + aut, y + dt * k2 * 0.5)
            k4 = dydt (t + dt,  y + dt * k3)
        
            ynew = y + dt * (k1 + (k2 + k3) * 2 + k4) * C1_6
        end subroutine Runge_Kutta4

        subroutine Gauss_Legendre4 (t, y, dt, dydt, ynew) ! Implicit
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y), 2)           :: rk

            rk = 0.
            call solve_rk_implicit (t, y, dt, dydt, FunK, rk)

            ynew = y + dt * (rk(:,1) + rk(:,2)) * 0.5
        
            contains
                
                subroutine FunK (t, y, dt, dydt, kin, kout)
                    implicit none
                    real*8, intent(in)                                          :: t, dt
                    procedure(dydt_tem)                                         :: dydt
                    real*8, dimension(:, :), intent(in)                         :: kin
                    real*8, dimension(size (kin, 1)), intent(in)                :: y
                    real*8, dimension(size (kin, 1), size (kin,2)), intent(out) :: kout

                    kout(:,1) = dydt (t + dt * (0.5 - SQ3_6), y + dt * ( &
                        & kin(:,1) * 0.25 + &
                        & kin(:,2) * (0.25 - SQ3_6)))
                    kout(:,2) = dydt (t + dt * (0.5 + SQ3_6), y + dt * ( &
                        & kin(:,1) * (0.25 + SQ3_6) + &
                        & kin(:,2) * 0.25))
                end subroutine FunK
        end subroutine Gauss_Legendre4 

        subroutine Runge_Kutta_four_oct4 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: k1, k2, k3, k4, kaux
              
            k1   = dydt (t,             y)
            kaux = k1 * C1_3
            k2   = dydt (t + dt * C1_3, y + dt * kaux)
            k3   = dydt (t + dt * C2_3, y + dt * (- kaux + k2))
            k4   = dydt (t + dt,        y + dt * (k1 - k2 + k3))
        
            ynew = y + dt * (k1 + (k2 + k3) * 3 + k4) * 0.125  
        end subroutine Runge_Kutta_four_oct4

        subroutine Runge_Kutta5 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: k1, k2, k3, k4, k5, k6
            real*8                                   :: aut
              
            aut = dt * 0.25
              
            k1 = dydt (t,             y)
            k2 = dydt (t + aut,       y + dt * k1 * 0.25)
            k3 = dydt (t + aut,       y + dt * (k1 + k2) * 0.125)
            k4 = dydt (t + dt * 0.5,  y + dt * (- k2 * 0.5 + k3))
            k5 = dydt (t + dt * 0.75, y + dt * (k1 + k4 * 3) * 3/16.)
            k6 = dydt (t + dt,        y + dt * (- k1 * 3 + k2 * 2 + (k3 - k4) * 12 + k5 * 8) * C1_7)
        
            ynew = y + dt * ((k1 + k6) * 7 + (k3 + k5) * 32 + k4 * 12)/90.
        end subroutine Runge_Kutta5

        subroutine Gauss_Legendre6 (t, y, dt, dydt, ynew) ! Implicit
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y), 3)           :: rk

            rk = 0.
            call solve_rk_implicit (t, y, dt, dydt, FunK, rk)

            ynew = y + dt * ((rk(:,1) + rk(:,3)) * 5 + rk(:,2) * 8) * C1_18
        
            contains
                
                subroutine FunK (t, y, dt, dydt, kin, kout)
                    implicit none
                    real*8, intent(in)                                          :: t, dt
                    procedure(dydt_tem)                                         :: dydt
                    real*8, dimension(:, :), intent(in)                         :: kin
                    real*8, dimension(size (kin, 1)), intent(in)                :: y
                    real*8, dimension(size (kin, 1), size (kin,2)), intent(out) :: kout

                    kout(:,1) = dydt (t + dt * (1. - SQ15_5) * 0.5,  y + dt * (&
                            & kin(:,1) * C5_36 + &
                            & kin(:,2) * (C2_3 - SQ15_5) * C1_3 + &
                            & kin(:,3) * (C5_36 - SQ15_30)))
                    kout(:,2) = dydt (t + dt * 0.5, y + dt * (&
                            & kin(:,1) * (C5_36 + SQ15_24)+ &
                            & kin(:,2) * C2_9 + &
                            & kin(:,3) * (C5_36 - SQ15_24)))
                    kout(:,3) = dydt (t + dt * (1. + SQ15_5) * 0.5, y + dt * (&
                            & kin(:,1) * (C5_36 + SQ15_30) + &
                            & kin(:,2) * (C2_9 + SQ15_15) + &
                            & kin(:,3) * C5_36))
                end subroutine FunK
        end subroutine Gauss_Legendre6 

        subroutine Runge_Kutta6 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: k1, k2, k3, k4, k5, k6, k7
              
            k1 = dydt (t,             y)
            k2 = dydt (t + dt * C1_3, y + dt * k1 * C1_3)
            k3 = dydt (t + dt * C2_3, y + dt * k2 * C2_3)
            k4 = dydt (t + dt * C1_3, y + dt * (k1 + k2 * 4 - k3) * C1_12)
            k5 = dydt (t + dt * C5_6, y + dt * (k1 * 25 - k2 * 110 + k3 * 35 + k4 * 90) * C1_48)
            k6 = dydt (t + dt * C1_6, y + dt * (k1 * 3 + k4 * 10 + k5 * 2) * 0.05 + (-k2 * 11 - k3 * 3)/24.)
            k7 = dydt (t + dt,        y + dt * (- k1 * 30.75 + k2 * 495 + k3 * 53.75 - k4 * 590 + k5 * 32 + k6 * 400)/195.)
        
            ynew = y + dt * ((k1 + k7) * 13/200. + (k3 + k4) * 11/40. + (k5 + k6) * 4/25.)
        end subroutine Runge_Kutta6

        subroutine Abbas6 (t, y, dt, dydt, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, dimension(size (y))              :: k1, k2, k3, k4, k5, k6, k7
            real*8                                   :: aut
              
            aut = dt * C1_3
              
            k1 = dydt (t,             y)
            k2 = dydt (t + aut,       y + dt * k1 * C1_3)
            k3 = dydt (t + dt * C2_3, y + dt * k2 * C2_3)
            k4 = dydt (t + aut,       y + dt * (k1 + k2 * 4 - k3) * C1_12)
            k5 = dydt (t + dt * C5_6, y + dt * (k1 * 25 - k2 * 110 + k3 * 35 + k4 * 90)/48.)
            k6 = dydt (t + dt * C1_6, y + dt * (k1 * 0.15 - k2 * 0.55 - k3 * 0.125 + k4 * 0.5 + k5 * 0.1))
            k7 = dydt (t + dt,        y + dt * (- k1 * 195.75 + k2 * 495 + k3 * 53.75 - k4 * 590 + k5 * 32 + k6 * 400)/195.)
        
            ynew = y + dt * ((k1 + k7) * 13 + (k3 + k4) * 55 + (k5 + k6) * 32) * 0.005
        end subroutine Abbas6

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
            
            ynew = y + dt * (k1 + k2) * 0.5
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
            k2 = dydt (t + dt * 0.5, y + dt * k1 * 0.5)

            ynew = y + dt * (k1 + k2 * 255) * 0.00390625

            k3 = dydt (t + dt, ynew)
            
            yaux = y + dt * (k1 + k2 * 500 + k3) * 0.001953125
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
            k2 = dydt (t + dt * 0.5,  y + dt * k1 * 0.5)
            k3 = dydt (t + dt * 0.75, y + dt * k2 * 0.75)

            ynew = y + dt * (k1 * 2 + k2 * 3 + k3 * 4) * C1_9

            k4 = dydt (t + dt, ynew)           
            
            yaux = y + dt * (k1 * 7/24. + k2 * 0.25 + k3 * C1_3 + k4 * 0.125)
        end subroutine Bogacki_Shampine3_2

        subroutine Zonneveld4_3 (t, y, dt, dydt, osol, oaux, yaux, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            integer*4, intent(out)                   :: osol, oaux
            real*8, dimension(size (y)), intent(out) :: ynew, yaux
            real*8, dimension(size (y))              :: k1, k2, k3, k4, k5, kaux
            real*8                                   :: aut
              
            aut = dt * 0.5
            
            osol = 4
            oaux = 3
            
            k1 = dydt (t,             y)
            k2 = dydt (t + aut,       y + dt * k1 * 0.5)
            k3 = dydt (t + aut,       y + dt * k2 * 0.5)
            k4 = dydt (t + dt,        y + dt * k3)
            k5 = dydt (t + dt * 0.75, y + dt * (k1 * 5 + k2 * 7 + k3 * 13 - k4) * 0.03125)
            
            kaux = (k2 + k3)

            ynew = y + dt * (k1 + kaux * 2 + k4) * C1_6
            yaux = y + dt * (- k1 * 3 + kaux * 14 + k4 * 13 - k5 * 32) * C1_6
        end subroutine Zonneveld4_3

        subroutine Merson4_3 (t, y, dt, dydt, osol, oaux, yaux, ynew)
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
            k2 = dydt (t + dt * C1_3, y + dt * k1 * C1_3)
            k3 = dydt (t + dt * C1_3, y + dt * (k1 + k2) * C1_6)
            k4 = dydt (t + dt * 0.5,  y + dt * (k1 + k3 * 3) * 0.125)
            k5 = dydt (t + dt,        y + dt * (k1 - k3 * 3 + k4 * 4) * 0.5)
            
            ynew = y + dt * (k1 + k4 * 4 + k5) * C1_6
            yaux = y + dt * (k1 + k3 * 3 + k4 * 4 + k5 * 2) * 0.1
        end subroutine Merson4_3

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
            k2 = dydt (t + dt * 0.25,   y + dt * k1 * 0.25)
            k3 = dydt (t + dt * 0.375,  y + dt * (k1 * 3 + k2 * 9) * 0.03125)
            k4 = dydt (t + dt * 12/13., y + dt * (k1 * 1932 - k2 * 7200 + k3 * 7296)/2197.)
            k5 = dydt (t + dt,          y + dt * ((k1 * 8341 + k3 * 29440 - k4 * 845)/4104. - k2 * 8))
            k6 = dydt (t + dt * 0.5,    y + dt * ((- k1 * 1216 + k4 * 1859)/4104. + k2 * 2 - k3 * 3544/2565. - k5 * 0.275))
            
            ynew = y + dt * ((k1 * 475 + k4 * 2197)/4104. + k3 * 1408/2565. - k5 * 0.2)
            yaux = y + dt * ((k1 * 6688 + k4 * 28561 + k6 * 2052)/56430. + k3 * 6656/12825. - k5 * 0.18)
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
            k2 = dydt (t + dt * 0.2,   y + dt * k1 * 0.2)
            k3 = dydt (t + dt * 0.3,   y + dt * (k1 * 3 + k2 * 9) * 0.025)
            k4 = dydt (t + dt * 0.6,   y + dt * (k1 * 3 - k2 * 9 + k3 * 12) * 0.1)
            k5 = dydt (t + dt,         y + dt * (- k1 * 11 + k2 * 135 - k3 * 140 + k4 * 70)/54.)
            k6 = dydt (t + dt * 0.875, y + dt * (k1 * 3262 + k2 * 37800 + k3 * 4600 + k4 * 44275 + k5 * 6831)/110592.)
            
            ynew = y + dt * (k1 * 37/378. + k3 * 250/621. + k4 * 125/594. + k6 * 512/1771.)
            yaux = y + dt * ((k1 * 5650 + k4 * 13525)/55296. + k3 * 18575/48384. + k5 * 277/14336. + k6 * 0.25)
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
            k2 = dydt (t + dt * 0.2,  y + dt * k1 * 0.2)
            k3 = dydt (t + dt * 0.3,  y + dt * (k1 * 0.075 + k2 * 0.225))
            k4 = dydt (t + dt * 0.8,  y + dt * (k1 * 44 - k2 * 168 + k3 * 160)/45.)
            k5 = dydt (t + dt * C8_9, y + dt * (k1 * 19372 - k2 * 76080 + k3 * 64448 - k4 * 1908)/6561.)
            k6 = dydt (t + dt,        y + dt * ((k1 * 9017 - k2 * 34080 + k4 * 882)/3168. + k3 * 46732/5247. - k5 * 5103/18656.))
            
            ynew = y + dt * ((k1 * 35 + k4 * 250)/384. + k3 * 500/1113. - k5 * 2187/6784. + k6 * 11/84.)
            
            k7 = dydt (t + dt, ynew)

            yaux = y + dt * ((k1 * 5179 + k4 * 35370)/57600. + k3 * 7571/16695. - k5 * 92097/339200. + k6 * 187/2100. + k7 * 0.025)
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
            k2 = dydt (t + dt * C1_6,  y + dt * k1 * C1_6)
            k3 = dydt (t + dt * C4_15, y + dt * (k1 * 4 + k2 * 16)/75.)
            k4 = dydt (t + dt * C2_3,  y + dt * (k1 * 5 - k2 * 16 + k3 * 15) * C1_6)
            k5 = dydt (t + dt * C5_6,  y + dt * ((- k1 * 165 - k3 * 425)/64. + (k2 * 880 + k4 * 85)/96.))
            k6 = dydt (t + dt,         y + dt * ((k1 * 612 + k5 * 88)/255. - k2 * 8 + (k3 * 4015 - k4 * 187)/612.))
            k7 = dydt (t + dt * C1_15, y + dt * ((- k1 * 8263 + k2 * 24800)/15000. - k3 * 643/680. - k4 * 0.324 + k5 * 2484/10625.))
            k8 = dydt (t + dt,         y + dt * (k1 * 3501/1720. + (297275 * k3 - 367200 * k2)/52632. - k4 * 319/2322. + &
                                         & k5 * 24068/84065. + k7 * 3850/26703.))
            
            ynew = y + dt * (k1 * 0.075 + k3 * 875/2244. + (k4 * 3703 + k7 * 125)/11592. + k5 * 264/1955. + k8 * 43/616.)
            yaux = y + dt * ((k1 * 13 + k4 * 50)/160. + (k3 * 2375 + k6 * 408)/5984. + k5 * 12/85.)
        end subroutine Verner6_5

        subroutine Fehlberg7_8 (t, y, dt, dydt, osol, oaux, yaux, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            integer*4, intent(out)                   :: osol, oaux
            real*8, dimension(size (y)), intent(out) :: ynew, yaux
            real*8, dimension(size (y))              :: k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, kaux
            
            osol = 7
            oaux = 8
            
            k1   = dydt (t,              y)
            k2   = dydt (t + dt * C2_27, y + dt * k1 * C2_27)
            k3   = dydt (t + dt * C1_9,  y + dt * (k1 + k2 * 3)/36.)
            k4   = dydt (t + dt * C1_6,  y + dt * (k1 + k3 * 3)/24.)
            k5   = dydt (t + dt * C5_12, y + dt * (k1 * C5_12 + (k4 - k3) * 1.5625))
            k6   = dydt (t + dt * 0.5,   y + dt * (k1 + k4 * 5 + k5 * 4) * 0.02)
            k7   = dydt (t + dt * C5_6,  y + dt * (- k1 * 25 + k4 * 125 - k5 * 260 + k6 * 250)/108.)
            k8   = dydt (t + dt * C1_6,  y + dt * (k1 * 93 + k5 * 244 - k6 * 200 + k7 * 13)/900.)
            k9   = dydt (t + dt * C2_3,  y + dt * (k1 * 2 + (- k4 * 795 + k5 * 1408 - k6 * 1070 + k7 * 67)/90. + k8 * 3))
            k10  = dydt (t + dt * C1_3,  y + dt * ((- k1 * 91 + k4 * 23  + k6 * 622)/108. - k5 * 976/135. + &
                                           & (- k7 * 19 + k8 * 170 - k9 * 5)/60.))
            kaux = - k4 * 8525 + k5 * 17984
            k11  = dydt (t + dt,         y + dt * (k1 * 2383 + kaux - k6 * 15050 + k7 * 2133 + &
                                           & k8 * 2250 + k9 * 1125 + k10 * 1800)/4100.)
            k12  = dydt (t,              y + dt * (((k1 - k7)/205. + ((- k6 + k10) * 2 + (- k8 + k9))/41.) * 3))
            k13  = dydt (t + dt,         y + dt * ((- k1 * 1777 + kaux - k6 * 14450 + k7 * 2193 + &
                                           & k8 * 2550 + k9 * 825 + k10 * 1900)/4100. + k12))

            kaux = k6 * 272 + (k7 + k8) * 216 + (k9 + k10) * 27

            ynew = y + dt * ((k1 + k11) * 41 + kaux) * C1_840
            yaux = y + dt * (kaux + (k12 + k13) * 41) * C1_840
        end subroutine Fehlberg7_8
        
        subroutine Dormand_Prince8_7 (t, y, dt, dydt, osol, oaux, yaux, ynew)
            implicit none
            real*8, intent(in)                       :: t, dt
            real*8, dimension(:), intent(in)         :: y
            procedure(dydt_tem)                      :: dydt
            integer*4, intent(out)                   :: osol, oaux
            real*8, dimension(size (y)), intent(out) :: ynew, yaux
            real*8, dimension(size (y))              :: k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13

            osol = 8
            oaux = 7

            k1  = dydt (t                               , y)
            k2  = dydt (t + dt *                   C1_18, y + dt * (                  C1_18 * k1))
            k3  = dydt (t + dt *                   C1_12, y + dt * (                  C1_48 * k1 + &
            & 0.0625 * k2))
            k4  = dydt (t + dt *                   0.125, y + dt * (                0.03125 * k1 + &
            &  0.09375 * k3))
            k5  = dydt (t + dt *                  0.3125, y + dt * (                 0.3125 * k1 - &
            & 1.171875 * k3 +  1.171875 * k4))
            k6  = dydt (t + dt *                   0.375, y + dt * (                 0.0375 * k1 + &
            &                    0.1875 * k4 +                     0.15 * k5))
            k7  = dydt (t + dt *                  0.1475, y + dt * (   29443841./614563906. * k1 + &
            &      77736538./692538347. * k4 -   28693883./1125000000. * k5 +    23124283./1800000000. * k6))
            k8  = dydt (t + dt *                   0.465, y + dt * (   16016141./946692911. * k1 + &
            &      61564180./158732637. * k4 +    22789713./633445777. * k5 +   545815736./2771057229. * k6 - &
            &   180193667./1043307555. * k7))
            k9  = dydt (t + dt * 5490023248./9719169821., y + dt * (   39632708./573591083. * k1 - &
            &     433636366./683701615. * k4 -  421739975./2616292301. * k5 +    100302831./723423059. * k6 + &
            &    790204164./839813087. * k7 +   800635310./3783071287. * k8))
            k10 = dydt (t + dt *                    0.65, y + dt * ( 246121993./1340847787. * k1 - &
            & 37695042795./15268766246. * k4 -  309121744./1061227803. * k5 -     12992083./490766935. * k6 + &
            &  6005943493./2108947869. * k7 +   393006217./1396673457. * k8 +   123872331./1001029789. * k9))
            k11 = dydt (t + dt * 1201146811./1299019798., y + dt * (-1028468189./846180014. * k1 + &
            &    8478235783./508512852. * k4 + 1311729495./1432422823. * k5 - 10304129995./1701304382. * k6 - &
            & 48777925059./3047939560. * k7 + 15336726248./1032824649. * k8 - 45442868181./3398467696. * k9 + &
            &  3065993473./597172653. * k10))
            k12 = dydt (t + dt                          , y + dt * (  185892177./718116043. * k1 - &
            &    3185094517./667107341. * k4 -  477755414./1098053517. * k5 -    703635378./230739211. * k6 + &
            &  5731566787./1027545527. * k7 +   5232866602./850066563. * k8 -   4093664535./808688257. * k9 + &
            & 3962137247./1805957418. * k10 +  65686358./487910083. * k11))
            k13 = dydt (t + dt                          , y + dt * (  403863854./491063109. * k1 - &
            &    5068492393./434740067. * k4 -   411421997./543043805. * k5 +    652783627./914296604. * k6 + &
            &  11173962825./925320556. * k7 - 13158990841./6184727034. * k8 +  3936647629./1978049680. * k9 - &
            &   160528059./685178525. * k10 + 248638103./1413531060. * k11))

            ynew  = y + dt * (13451932./455176623. * k1 - 808719846./976000145. * k6 + &
            & 1757004468./5645159321. * k7 + 656045339./265891186. * k8 - 3867574721./1518517206. * k9 + &
            &  465885868./322736535. * k10 +   53011238./667516719. * k11 +                  C2_45 * k12)
            yaux  = y + dt * (14005451./335480064. * k1 - 59238493./1068277825. * k6 + &
            &   181606767./758867731. * k7 + 561292985./797845732. * k8 - 1041891430./1371343529. * k9 + &
            & 760417239./1151165299. * k10 +  118820643./751138087. * k11 - 528747749./2220607170. * k12 + 0.25 * k13)
        end subroutine Dormand_Prince8_7

        !!

        !! Runge Kutta Half_Step

        recursive subroutine rk_half_step (t, y, dt_adap, dydt, integ, ord, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            integer*4, intent(in)                    :: ord
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            procedure(integ_tem)                     :: integ
            real*8, dimension(:), intent(in)         :: y
            real*8, dimension(size (y))              :: yhalf, yaux, yscal
            procedure(dydt_tem)                      :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, intent(inout)                    :: dt_adap, dt_used
            real*8                                   :: e_calc, ratio, hdt_adap
            integer*4                                :: iter = 0

            iter     = iter + 1
            dt_adap  = max (dt_adap, dt_min)
            hdt_adap = 0.5 * dt_adap

            call integ (           t,     y,  dt_adap, dydt,  ynew)
            call integ (           t,     y, hdt_adap, dydt, yhalf)                       
            call integ (t + hdt_adap, yhalf, hdt_adap, dydt,  yaux)

            yscal = abs (y + dt_adap * dydt (t, y)) + SAFE_LOW
            
            e_calc = max (maxval (abs ((ynew - yaux) / yscal) / real(2**ord - 1, kind=8)), SAFE_LOW)
            ratio  = e_tol / e_calc
            if (ratio > 1.) then
                dt_used = dt_adap
                dt_adap = dt_adap * min (beta * ratio**(1. / real (ord + 1, kind=8)), MAX_DT_FAC)
                iter    = 0
            else
                if (dt_adap .eq. dt_min) then
                    dt_used = dt_min
                    iter = 0
                else
                    dt_adap = dt_adap * min (beta * ratio**(1. / real (ord, kind=8)), MAX_DT_FAC)
                    if ((isnan (dt_adap)) .or. (dt_adap .le. dt_min) .or. (iter .eq. MAX_N_ITER)) then
                        dt_used = dt_min
                        dt_adap = dt_min
                        call integ (t, y, dt_adap, dydt, ynew)
                        iter = 0
                    else
                        call rk_half_step (t, y, dt_adap, dydt, integ, ord, e_tol, beta, dt_min, dt_used, ynew)
                    end if
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
        subroutine rk_adap_wrapper (t, y, dt_adap, dydt, integ, ord, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            integer*4, intent(in)                    :: ord
            real*8, intent(in)                       :: t, e_tol, beta, dt_min
            procedure(integ_tem)                     :: integ
            real*8, dimension(:), intent(in)         :: y
            type(dydt_i), dimension(size (y))        :: dydt
            real*8, dimension(size (y)), intent(out) :: ynew
            real*8, intent(inout)                    :: dt_adap, dt_used
            
            call rk_half_step_caller (t, y, dt_adap, Faux, integ, ord, e_tol, beta, dt_min, dt_used, ynew)

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
