! space_time_predictor.f90
! author: sunder

!-----------------------------------------------------------------------
! Apply ADER-DG method to the reconstruction
!-----------------------------------------------------------------------

subroutine ader_space_time_predictor(lqhi, lqbnd, lFbnd, luh)

    use ader_weno
    implicit none
    ! Argument list
    double precision, intent(in)  :: luh(nVar,nDOF)                ! spatial degrees of freedom
    double precision, intent(out) :: lqhi(nVar,nDOF)               ! time-averaged space-time degrees of freedom
    double precision, intent(out) :: lqbnd(nVar,2), lFbnd(nVar,2)  ! time-averaged space-time degrees of freedom
    ! Local variables
    integer :: i,l,iVar, iter
    double precision    :: lFhi(nVar,nDOF)            ! time-averaged nonlinear flux tensor in each space-time DOF
    double precision    :: rhs0(nVar,nDOF,nDOF)       ! contribution of the initial condition to the known right hand side
    double precision    :: rhs(nVar,nDOF,nDOF)        ! known right hand side
    double precision    :: lqh(nVar,nDOF, nDOF)       ! space-time degrees of freedom
    double precision    :: lFh(nVar,nDOF,nDOF)        ! nonlinear flux tensor in each space-time DOF
    double precision    :: lqhold(nVar,nDOF,nDOF)     ! old space-time degrees of freedom
    
    do i = 1, nDOF
        ! Trivial initial guess (can be significantly improved)
        do iVar = 1, nVar
            lqh(iVar,i,:) = luh(iVar,i)
        end do
        ! Compute the contribution of the initial condition uh to the time update.
        do iVar = 1, nVar
            rhs0(iVar,i,:) = wGP(i)*F0(:)*luh(iVar,i)
        end do
    end do
    
    ! Discrete Picard iterations.
    
    do iter = 1, N+1
        
        lqhold = lqh ! Save old space-time DOF
        
        do l = 1, nDOF ! loop over DOF in time
            
            ! Compute the fluxes (once these fluxes are available, the subsequent operations are independent from each other)
            
            do i = 1, nDOF
                call PDEFlux(lFh(:,i,l),lqh(:,i,l))
            end do
            
            ! Compute the "derivatives" (contributions of the stiffness matrix)

            rhs(:,:,l) = rhs0(:,:,l) - wGP(l)*(dt/dx)*MATMUL( lFh(:,:,l), Kxi )
        
        end do ! end loop over time DOF
        
        
        ! Multiply with (K1)^(-1) to get the discrete time integral of the discrete Picard iteration
        
        do i = 1, nDOF
             lqh(:,i,:) = 1.0d0/(wGP(i))*MATMUL( rhs(:,i,:), TRANSPOSE(iK1) )
        end do
        
    end do
    
    ! Immediately compute the time-averaged space-time polynomials
    
    do i = 1, nDOF
        lqhi(:,i) = matmul( lqh(:,i,:), wGP )
        lFhi(:,i) = matmul( lFh(:,i,:), wGP )
    end do

    ! Compute the bounday-extrapolated values for Q and F
    
    lqbnd(:,1) = matmul( lqhi(:,:), phiL )   ! left Q
    lqbnd(:,2) = matmul( lqhi(:,:), phiR )   ! right Q
    lFbnd(:,1) = matmul( lFhi(:,:), phiL )   ! left F
    lFbnd(:,2) = matmul( lFhi(:,:), phiR )   ! right F

end subroutine ader_space_time_predictor
