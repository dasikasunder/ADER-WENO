! write_data.f90
! author: sunder

!-----------------------------------------------------------------------
! Output data as csv file
!-----------------------------------------------------------------------

subroutine write_data
    use ader_weno
    implicit none
    integer :: i
    
    open(7,file='sol.csv')
    
    write(7,'(a,a,a)') 'x', ',', 'U'
    do i = 1, IMAX
        write(7,'(es14.7,a,es14.7)') x(i), ',', uh(1,i)
    end do
    
    close(7)

end subroutine write_data
