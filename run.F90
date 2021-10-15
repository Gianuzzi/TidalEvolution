module run
    
    implicit none
    integer*8     :: n_iter                 ! N_iterations
    integer       :: n_points               ! N_outputs
    integer*8     :: i                      ! LOOP
    real*8        :: t, t0, tf, dt, dt_adap ! Times
    real*8        :: Logt, t_add, t_out     ! Time [output]
    character(32) :: filename               ! Output file
    character(1)  :: selection              ! In case the file exists
    logical       :: file_exists            ! File existance  
    real*8        :: start_time             ! Excecution start time
    real*8        :: final_time             ! Total elapsed time

end module run
