! read_input.f90
! author: sunder

!-----------------------------------------------------------------------
! Read the input data from the parameters file
!-----------------------------------------------------------------------

subroutine read_input
    use ader_weno
    implicit none
    
    open(UNIT = 1, FILE = 'params.ini', STATUS = 'UNKNOWN')
    
    read(1,*)EQN%a
    read(1,*)N
    read(1,*)IMAX
    read(1,*)tend
    read(1,*)xL
    read(1,*)xR
    read(1,*)bL
    read(1,*)bR
    read(1,*)ICType
    
    close(1)
    
end subroutine
