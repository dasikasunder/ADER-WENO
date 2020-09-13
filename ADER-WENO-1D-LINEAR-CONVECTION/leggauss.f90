! leggaus.f90
! author: sunder

! ---------------------------------------------------------------------------------------------------------
! Gauss-Legendre quadrature points (x) and weights (w) in interval (x1, x2) for n points   
! ---------------------------------------------------------------------------------------------------------

pure subroutine gauleg(x1, x2, x, w, n)
    use constants, only : m_pi
    implicit none
    ! Argument list
    integer, intent(in)           ::  n
    double precision, intent (in) :: x1, x2
    double precision, intent(out) :: x(n), w(n)
    ! Local variables
    double precision :: EPS
    integer  :: i,j,m
    double precision :: p1, p2, p3, pp, xl, xm, z, z1

    parameter (EPS=1.0D-15)

    m  = (n+1)/2
    xm = 0.5*(x2+x1)
    xl = 0.5*(x2-x1)
    do i=1,m
        z = COS(m_pi*(i-0.25d0)/(n+0.5d0))
    1    continue
        p1 = 1.0d0
        p2 = 0.0d0
        do j = 1,n
            p3 = p2
            p2 = p1
            p1 = ((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/dble(j)
        end do
        pp = dble(n)*(z*p1-p2)/(z*z-1.0d0)
        z1 = z
        z  = z1-p1/pp
        if(abs(z-z1).GT.EPS)goto 1
        x(i)    = xm-xl*z
        x(n+1-i)= xm+xl*z
        w(i)    = 2.0d0*xl/((1.-z*z)*pp*pp)
        w(n+1-i)= w(i)
    end do
    return
end subroutine gauleg
