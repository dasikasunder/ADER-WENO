! analyze_soln.f90
! author: sunder

!-----------------------------------------------------------------------
! Find the LMAX error of the solution for peridic case and print it
!-----------------------------------------------------------------------

subroutine analyze_soln
    use ader_weno
    implicit none
    integer :: i, q
    double precision :: u_exact, xp
    double precision :: q0(nVar)
    double precision :: error(IMAX)
    
    do i = 1, IMAX
        u_exact = 0.0 
        do q = 1, N+1
            xp = x(i) - 0.5*dx + dx*xGP(q)
            call initial_field(q0, xp)
            u_exact = u_exact + wGP(q)*q0(1)
        end do
        
        error(i) = abs(u_exact - uh(1,i))
        
    end do
    
    print *, 'Lmax = ', maxval(error)

end subroutine analyze_soln
