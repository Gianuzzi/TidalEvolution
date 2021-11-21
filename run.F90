module run
    
    implicit none
    integer*8                         :: n_iter                         ! N_iterations
    integer*8                         :: n_points                       ! N_outputs
    integer*8                         :: i, j, l                        ! LOOPs
    real*8                            :: t, t0, tf                      ! Times
    real*8                            :: dt, dt_adap, dt_min            ! dTimes
    real*8                            :: t_add, logt, fixt              ! Times [output]
    real*8, dimension(:), allocatable :: t_out                          ! Times [output]
    character(32)                     :: filename                       ! Output file
    character(1)                      :: selection                      ! In case the file exists
    logical, parameter                :: TRUE = .true., FALSE = .false. !True | False
    logical                           :: file_exists, logsp             ! File existance / LogT  
    real*8                            :: start_time, final_time         ! Excecution start and elapsed time
    real*8, parameter                 :: MIN_VAL = 1e-15                ! Avoid extremely low values
    real*8                            :: aux1, aux2, aux3               ! Aux variables

    contains

        subroutine get_t_outs(t0, tf, n_points, logsp, t_out)
            implicit none
            real*8, intent(in)                 :: t0, tf
            real*8, dimension(:), allocatable  :: t_out
            integer*8, intent(in)              :: n_points
            logical                            :: logsp
            real*8                             :: t, t_aux, t_add
            integer*8                          :: i

            if (n_points < 2) then
                write(*,*) "n_points should be greater than 2."
                stop (1)
            end if
            allocate (t_out(n_points))
            t_out(1) = t0
            t        = t0
            if (logsp .eqv. .TRUE.) then
                t_aux = exp (log (tf - t0 + 1.) / (real (n_points, kind=8)))
                t_add = t_aux
                do i = 2, n_points - 1
                    t_add    = t_add * t_aux
                    t_out(i) = t0 + t_add - 1.
                end do
            else
                t_aux = (tf - t0) / (real (n_points, kind=8) - 1.)
                do i = 1, n_points - 1
                    t_out(i + 1) = t0 + t_aux * i
                end do
            end if
            t_out(n_points) = tf
        end subroutine get_t_outs

end module run
