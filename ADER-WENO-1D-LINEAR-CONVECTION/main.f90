! main.f90
! author: sunder

!-----------------------------------------------------------------------
! Main program of the code. Puts together everything.
!-----------------------------------------------------------------------

program main
    use ader_weno
    implicit none

    
    call read_input
    call initialize
    
    do while (time < tend)
    
        if ((time + dt) .gt. tend) then
            dt = tend - time;
        end if
        
        call compute_rhs
        
        uh(:, 1:IMAX) = uh(:, 1:IMAX) + dt*duh(:,:)
        
        print *, 'dt = ', dt, 'time = ', time
        
        time = time + dt
        timestep = timestep + 1
        
    end do
    
    print *, 'dt = ', dt, 'time = ', time

    call analyze_soln
    call write_data
    
    print *, 'Done. Number of timesteps = ', timestep
    
end program main

