! compute_rhs.f90
! author: sunder

!-----------------------------------------------------------------------
! For a given stencil, find the oscillation indicator
!-----------------------------------------------------------------------

subroutine os_indicator(ul, OS)
    use ader_weno, only : N, OS_M
    implicit none
    ! Argument list 
    double precision, intent(in) :: ul(N+1)
    double precision, intent(out) :: OS
    ! Local variables
    integer :: p, m
    
    OS = 0.0d0
    
    do p = 1, N+1
        do m = 1, N+1
            OS = OS + OS_M(p,m)*ul(p)*ul(m)
        end do
    end do
    
end subroutine os_indicator

!-----------------------------------------------------------------------
! Mapping function of Henrick, Aslam and Powers
!-----------------------------------------------------------------------

function OS_map(non_lin_wt, lin_wt)

    implicit none
    double precision :: OS_map
    double precision, intent(in) :: non_lin_wt, lin_wt
    double precision :: num, denom

    num = non_lin_wt*(lin_wt + lin_wt*lin_wt - 3.0d0*lin_wt*non_lin_wt + non_lin_wt*non_lin_wt)
    denom = lin_wt*lin_wt + non_lin_wt*(1.0d0-2.0d0*lin_wt)

    OS_map = num/denom

end function OS_map

!-----------------------------------------------------------------------
! Given stencil in 1D, reconstruct the solution
!-----------------------------------------------------------------------

subroutine reconstruct_1d(lu, luh)
    use ader_weno
    implicit none
    ! Argument list
    double precision, intent(in) :: lu(-N:N)  ! Stencil containing average values
    double precision, intent(out) :: luh(N+1) ! Unknown coefficients to be found
    ! Local variables
    integer :: L, iDOF, i
    double precision :: ustencil(N+1)         ! Local stencil
    double precision, dimension(N+1) :: luh_l, luh_r, luh_cl, luh_cr
    double precision :: OS_L, OS_R, OS_CL, OS_CR
    double precision :: w(nStencils), total, OS_map
    logical, parameter :: map = .false.

    ! Left biased stencil

    L = -N;

    ustencil = lu(L:L+N)
    luh_l = matmul(iML, ustencil)
    call os_indicator(luh_l, OS_L)
    w(1) = lin_wt(1)/(OS_L + small_num)**4

    ! Right biased stencil

    L = 0

    ustencil = lu(L:L+N)
    luh_r = matmul(iMR, ustencil)
    call os_indicator(luh_r, OS_R)
    w(2) = lin_wt(2)/(OS_R + small_num)**4

    ! Center right stencil

    L = -floor(real(N/2))

    ustencil = lu(L:L+N)
    luh_cr = matmul(iMCR, ustencil)
    call os_indicator(luh_cr, OS_CR)
    w(3) = lin_wt(3)/(OS_CR + small_num)**4

    if (nStencils .eq. 4) then
        L = -(floor(real(N/2)) + 1)

        ustencil = lu(L:L+N)
        luh_cl = matmul(iMCL, ustencil)
        call os_indicator(luh_cl, OS_CL)
        w(4) = lin_wt(4)/(OS_CL + small_num)**4

    end if

    ! Normalize the weights

    total = sum(w)
    w = w/total

    if (map) then

        ! Find the maped weights

        do i = 1, nStencils
            w(i) = OS_map(w(i), lin_wt(i))
        end do

        total = sum(w)
        w = w/total

    end if

    do idoF = 1, N+1
        luh(iDOF) = w(1)*luh_l(iDOF) + w(2)*luh_r(iDOF) + w(3)*luh_cr(iDOF)
        if(nStencils .eq. 4) then
            luh(iDOF) = luh(iDOF) + w(4)*luh_cl(iDOF)
        end if
    end do

end subroutine reconstruct_1d

!-----------------------------------------------------------------------
! Find the update coefficients for each cell
!-----------------------------------------------------------------------

subroutine compute_rhs
    use ader_weno
    implicit none

    double precision :: lu(-N:N), luh(n+1), uhi(nVar,N+1), qhi(nVar,N+1)
    integer :: i, k, iVar

    call apply_boundary_conditions

    ! Find boundary extrapolated values of conserved variables in each cell

    do i = 0, IMAX+1
        do iVar = 1, nVar

            do k = -N, N
                lu(k) = uh(iVar, i+k)
            end do

            call reconstruct_1d(lu, luh)

            do k = 1, N+1
                uhi(iVar,k) = luh(k)
            end do

        end do

        call ader_space_time_predictor(qhi, qbnd(:,:,i), Fbnd(:,:,i), uhi)

    end do

    ! Solve riemann problem on each face

    do i = 1, IMAX + 1
        call RiemannFlux(F(:,i), qbnd(:,2,i-1), qbnd(:,1,i), Fbnd(:,2,i-1), Fbnd(:,1,i))
    end do

    ! Find the update values in each cell

    do i = 1, IMAX
        duh(:,i) = -(1.0d0/dx)*( F(:,i+1) - F(:,i) )
    end do

end subroutine compute_rhs
