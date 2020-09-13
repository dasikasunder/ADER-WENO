! pde.f90
! author: sunder

!-----------------------------------------------------------------------
! Convert conserved variable to primitive variable
!-----------------------------------------------------------------------

subroutine Cons2Prim(V, Q)
    use ader_weno, only : nVar
    implicit none
    
    ! Argument list 
    double precision, intent(in)   :: Q(nVar)
    double precision, intent(out)  :: V(nVar)

    V(1) = Q(1)

end subroutine Cons2Prim

!-----------------------------------------------------------------------
! Convert primitive variable to conserved variable
!-----------------------------------------------------------------------

subroutine Prim2Cons(Q, V)
    USE ader_weno, ONLY : nVar
    implicit none
    
    ! Argument list 
    double precision, intent(IN)   :: V(nVar)
    double precision, intent(OUT)  :: Q(nVar)
    
    Q(1) = V(1)
    
end subroutine Prim2Cons

!-----------------------------------------------------------------------
! Find the conservative flux for the PDE
!-----------------------------------------------------------------------

subroutine PDEFlux(F,Q)
    USE ader_weno, ONLY : nVar, EQN
    implicit none
    
    ! Argument list 
    double precision, intent(IN)   :: Q(nVar)
    double precision, intent(OUT)  :: F(nVar)
    
    F(1) = EQN%a*Q(1)

end subroutine PDEFlux

!-----------------------------------------------------------------------
! Find the eigenvalues of the PDE
!-----------------------------------------------------------------------

subroutine PDEEigenvalues(Lambda, Q)
    use ader_weno, only : nVar, EQN
    implicit none
    ! Argument list
    double precision, intent(IN)  :: Q(nVar)
    double precision, intent(OUT) :: Lambda(nVar)
    
    Lambda = (/ EQN%a /)                          ! Eigenvalues of the linear convection equations
    
end subroutine PDEEigenvalues

!-----------------------------------------------------------------------
! Compute the upwind riemann flux
!-----------------------------------------------------------------------

subroutine RiemannFlux(F, QL, QR, FL, FR)
    
    use ader_weno, ONLY : nVar
    implicit none
    
    ! Argument list declaration
    double precision, intent(IN)  :: QL(nVar), QR(nVar), FL(nVar), FR(nVar)
    double precision, intent(OUT) :: F(nVar)
    ! Local variable declaration
    double precision :: smax, LL(nVar), LR(nVar)

    call PDEEigenvalues(LL,QL)
    call PDEEigenvalues(LR,QR)
    
    smax = max( maxval(abs(LL)), maxval(abs(LR)) )
    
    F = 0.5d0*(FR + FL - smax*(QR - QL))
 
end subroutine RiemannFlux
