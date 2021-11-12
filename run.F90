module run
    
    implicit none
    integer*8         :: n_iter               ! N_iterations
    integer*8         :: n_points             ! N_outputs
    integer*8         :: i, j, l              ! LOOPs
    real*8            :: t, t0, tf            ! Times
    real*8            :: dt, dt_adap, dt_min  ! dTimes
    real*8            :: Logt, t_add, t_out   ! Times [output]
    character(32)     :: filename             ! Output file
    character(1)      :: selection            ! In case the file exists
    logical           :: file_exists          ! File existance  
    real*8            :: start_time           ! Excecution start time
    real*8            :: final_time           ! Total elapsed time
    real*8, parameter :: MIN_VAL = 1e-15      ! Avoid extremely low values
    real*8            :: aux1, aux2, aux3     ! Aux variables

end module run
