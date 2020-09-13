! initialize.f90
! author: sunder

!-----------------------------------------------------------------------
! Allocate memory, initialize various data and also the initial
! condition
!-----------------------------------------------------------------------

subroutine initialize
    use ader_weno
    implicit none
    
    integer :: i, j, p, m, alpha, q, center, L, iGP
    double precision, dimension(0:N) :: phi_p, phi_m, phi
    double precision, dimension(N+1,N+1) :: MS, FRm
    double precision, dimension(nVar) :: q0
    double precision, dimension(N+1) :: phi_xi, phi_i
    double precision :: xiL, xiR, shift, xp, total
    
    ! Initialize basic stuff 
    
    dx = (xR - xL)/REAL(IMAX)
    dt = CFL*dx 
    time = 0.0d0
    timestep = 0
    nDOF = N+1
    
    if (mod(N,2) .eq. 0) then
        nStencils = 3 
    else
        nStencils = 4
    end if
    
    nGhostCells = 2*(N+1) -1
    
    ! Allocate required memory 
    
    allocate(  x(1-nGhostCells:IMAX+nGhostCells) )
    allocate(  uh(nVar, 1-nGhostCells:IMAX+nGhostCells) )
    allocate(  duh(nVar, IMAX) )
    allocate(  qbnd(nVar, 2, 0:IMAX+1) ) ! 2 corresponds to 2 faces per cell
    allocate(  Fbnd(nVar, 2, 0:IMAX+1) )
    allocate(  F(nVar, IMAX+1) )
    
    ! Initialize the grid 
    
    do i = 1-nGhostCells, IMAX+nGhostCells
        x(i) = xL + (REAL(i)-0.5)*dx
    end do
    
    ! Basis functions 
    
    xiL = 0.0d0
    xiR = 1.0d0
    
    allocate(  xGP(N+1) )
    allocate(  wGP(N+1) )
    allocate(  phiL(N+1) )
    allocate(  phiR(N+1) )
    
    call gauleg(xiL, xiR, xGP, wGP, N+1)
    
    ! Find the basis function values at the left and the right end of the unit cell 
    
    do i = 1, N+1
        call lagrange_basis_values(i, xiL, phi)
        phiL(i) = phi(0)
        call lagrange_basis_values(i, xiR, phi)
        phiR(i) = phi(0)
    end do
    
    ! WENO Related data 
    
    allocate (lin_wt(nStencils)) ! Calculate linear weights

    lin_wt(1) = lambda_s
    lin_wt(2) = lambda_s
    lin_wt(3) = lambda

    if (nStencils .eq. 4) then
        lin_wt(4) = lambda
    end if

    total = sum(lin_wt)

    lin_wt = lin_wt/total

    allocate( OS_M(N+1, N+1) ) ! Initialize the oscillation indicator matrix
    
    OS_M = 0.0d0
    
    do p = 1, N+1
        do m = 1, N+1
            do q = 1, N+1
                call lagrange_basis_values(p, xGP(q), phi_p)
                call lagrange_basis_values(m, xGP(q), phi_m)
                do alpha = 1, N
                    OS_M(p,m) = OS_M(p,m) + wGP(q)*phi_p(alpha)*phi_m(alpha)
                end do
            end do
        end do
    end do
    
    allocate( iMCL(N+1, N+1) ) ! Initialize center left stencil coefficient matrix
    allocate( iML(N+1, N+1) )  ! Initialize left-sided coefficient matrix
    allocate( iMCR(N+1, N+1) ) ! Initialize center right-sided coefficient matrix
    allocate( iMR(N+1, N+1) )  ! Initialize right-sided coefficient matrix
    
    ! Centered left side stencil 
    
    MS = 0.0d0
    L = floor(real(N/2)) + 1
    
    do i = 1, N+1
        center = i - 1 - L
        do j = 1, N+1
            do q = 1, N+1
                shift = xGP(q) + DBLE(center)
                call lagrange_basis_values(j, shift, phi)
                MS(i,j) = MS(i,j) + wGP(q)*phi(0)
            end do
        end do
    end do
    
    call MatrixInverse(N+1,MS,iMCL)
    
    ! Centered right side stencil (Use this stencil for odd-ordered schemes)
    
    MS = 0.0d0
    L = floor(real(N/2))
    
    do i = 1, N+1
        center = i - 1 - L
        do j = 1, N+1
            do q = 1, N+1
                shift = xGP(q) + DBLE(center)
                call lagrange_basis_values(j, shift, phi)
                MS(i,j) = MS(i,j) + wGP(q)*phi(0)
            end do
        end do
    end do
    
    call MatrixInverse(N+1,MS,iMCR)
    
    ! Left sided stencil 
    
    MS = 0.0d0
    L = N 
    
    do i = 1, N+1
        center = i - 1 - L
        do j = 1, N+1
            do q = 1, N+1
                shift = xGP(q) + dble(center)
                call lagrange_basis_values(j, shift, phi)
                MS(i,j) = MS(i,j) + wGP(q)*phi(0)
            end do
        end do
    end do
    
    call MatrixInverse(N+1,MS,iML)
    
    ! Right sided stencil 
    
    MS = 0.0d0
    L = 0
    
    do i = 1, N+1
        center = i - 1 - L
        do j = 1, N+1
            do q = 1, N+1
                shift = xGP(q) + DBLE(center)
                call lagrange_basis_values(j, shift, phi)
                MS(i,j) = MS(i,j) + wGP(q)*phi(0)
            end do
        end do
    end do
    
    call MatrixInverse(N+1,MS,iMR)
    
    ! ADER related data 
    
    allocate( Kxi(N+1, N+1) ) ! Element stiffness matrix
    allocate( K1(N+1, N+1) ) ! Element stiffness matrix
    allocate( iK1(N+1, N+1) ) ! Element stiffness matrix
    allocate( F0(N+1) ) ! Time flux matrix
    
    Kxi = 0.0d0
    
    DO iGP = 1, N+1
        call BaseFunc1D(phi_i,phi_xi,xGP(iGP))
        DO i = 1, N+1
            DO j = 1, N+1
                Kxi(i,j) = Kxi(i,j) + wGP(iGP)*phi_xi(i)*phi_i(j)
            end do
        end do
    end do
    
    do i = 1, N+1
        do j = 1, N+1
            FRm(i,j) = phiR(i)*phiR(j)   ! Left contribution to the right flux matrix   (m = left  of the interface)
        end do
    end do
    
    ! The time flux matrices for the ADER-DG predictor method are given by the principle of upwinding in time (causality principle)
    K1 = FRm - Kxi ! upwinding in time = information comes from smaller times
    F0 = phiL   ! upwinding in time = information comes from smaller times
    call MatrixInverse(N+1,K1,iK1)
    
    
    ! Initialize the solution 
    
    do i = 1, IMAX
        uh(:,i) = 0.0
        do q = 1, N+1
            xp = x(i) - 0.5d0*dx + dx*xGP(q)
            call initial_field(q0, xp)
            uh(:,i) = uh(:,i) + wGP(q)*q0
        end do
    end do
    
end subroutine initialize
