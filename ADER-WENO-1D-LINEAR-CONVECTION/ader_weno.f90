! ader_weno.f90
! author: sunder

!-----------------------------------------------------------------------
! Main module containing all the necessary data structures.
!-----------------------------------------------------------------------

module ader_weno
    use constants, only : m_pi
    implicit none
    public
    
    ! ---------------------------------------------------------------------------------------------------------
    ! Input parameters - This part can be modified by the user  
    ! ---------------------------------------------------------------------------------------------------------
    
    integer, parameter  :: nVar   = 1          ! Number of variables in the PDE system
    integer             :: N                   ! Degree of approximation
    integer             :: IMAX                ! Number of cells in the x-direction
    double precision    :: tend                ! End time of the simulation
    double precision    :: xL                  ! Left end of the domain
    double precision    :: xR                  ! Right end of the domain
    integer             :: bL                  ! Boundary condition on left  : 0-I, 1-T, 2-R, 3-P
    integer             :: bR                  ! Boundary condition on right : 0-I, 1-T, 2-R, 3-P
    double precision, parameter :: CFL = 0.9d0 ! CFL constant (should be less than 1.0)
    integer             :: ICType              ! Initial condition function.
                                               ! Optional -> Can be set in the initial field file as well in case desired test case does not exist
    
    ! Supported boundary conditions: I-inflow, T-transmissive, R-reflective, P-periodic 
    
    ! ---------------------------------------------------------------------------------------------------------
    ! Do not change anything below   
    ! ---------------------------------------------------------------------------------------------------------
    
    ! Basic stuff 
    
    double precision :: dx                           ! Size of each cell in the domain (assumed constant)
    double precision :: dt                           ! Time step size
    double precision :: time                         ! Current time in the simulation
    integer          :: timestep                     ! Current time step number
    integer          :: nDOF                         ! Number of degrees of freedom in each cell
    integer          :: nStencils                    ! Number of stencils for WENO reconstruction
    integer          :: nGhostCells                  ! Number of ghost cells on each side of the domain to apply boundary conditions
    
    ! Arrays holding the solution, grid, rhs and fluxes 
    
    double precision, allocatable  :: x(:)           ! Mesh - contains the cell center values
    double precision, allocatable  :: uh(:,:)        ! Solution vector (nVar, nCells+2*nGhostCells) -> also contains ghost values
    double precision, allocatable  :: duh(:,:)       ! RHS for updating the solution
    double precision, allocatable  :: qbnd(:,:,:)    ! Conserved variable on cell boundaries
    double precision, allocatable  :: Fbnd(:,:,:)    ! Conserved fluxes on cell boundaries
    double precision, allocatable  :: F(:,:)         ! Upwind flux at each face
    
    ! Basis functions 
    
    double precision, allocatable  :: xGP(:), wGP(:)   ! Legendre-Gauss quadrature points in interval [0,1]
    double precision, allocatable  :: phiL(:), phiR(:) ! Basis function values on the left and side of unit cell
    
    ! WENO related data 
    
    double precision, parameter   :: lambda_s  = 1.0d0   ! Weight for the sided stencils
    double precision, parameter   :: lambda    = 100.0d0 ! Weight for the centered stencils
    double precision, allocatable :: lin_wt(:)           ! Normalized linear weights of each stencil
    double precision, parameter   :: small_num = 1.0d-12 ! Small number (to avoid division by zero)
    double precision, allocatable :: OS_M(:,:)           ! Oscillation indicator matrix
    double precision, allocatable :: iML(:,:), iMCL(:,:) ! Inverse coefficient matrix for left and center left stencil
    double precision, allocatable :: iMR(:,:), iMCR(:,:) ! Inverse coefficient matrix for right and center right stencil
    
    ! ADER related data 
    double precision, allocatable :: Kxi(:,:)            ! Element stiffness matrix
    double precision, allocatable :: K1(:,:), iK1(:,:)   ! F1 - Ktau and its inverse
    double precision, allocatable :: F0(:)               ! Time flux matrices
    
    ! Important info and parameters concerning the governing PDE system

    type tEquations
        double precision :: a ! Advection speed of the wave
    end type tEquations

    type(tEquations)   :: EQN

end module ader_weno
