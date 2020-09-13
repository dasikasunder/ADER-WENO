! Lagrange polynomial basis 

subroutine lagrange_basis_values(center, xi, phi)
    use ader_weno
    implicit none
    ! Argument list
    integer, intent(IN) :: center
    double precision, intent(IN) :: xi
    double precision, intent(OUT) :: phi(0:N)
    ! Local variables 
    integer :: i, k
    double precision :: lagrange_weight, tmp_lagrange_weight, v, k_faculty
    
    tmp_lagrange_weight = 1.0d0
    k_faculty = 1.0d0

    ! First find the denominator 
    
    do i = 1, N+1
        if (i .ne. center) then
            tmp_lagrange_weight = tmp_lagrange_weight*(xGP(center) - xGP(i))
        end if
    end do
    
    lagrange_weight = 1.0d0/tmp_lagrange_weight
    
    phi(0) = 1.0d0
    
    do i = 1, N
        phi(i) = 0.0d0
    end do
    
    do i = 0, N
        if (i+1 .ne. center) then
            v = xi - xGP(i+1)
            do k = N, 1, -1
                phi(k) = (phi(k)*v + phi(k - 1));
            end do
            phi(0) = phi(0)*v;
        end if
    end do
    
    do k = 0, N
        phi(k) = phi(k)*k_faculty*lagrange_weight
        k_faculty = k_faculty*dble(k+1)
    end do
    
end subroutine lagrange_basis_values


subroutine BaseFunc1D(phi,phi_xi,xi)
   use ader_weno
   implicit none
   ! Argument list
   double precision, intent(in ) :: xi                              ! coordinate in [0,1] where to evaluate the basis
   double precision, intent(out) :: phi(N+1), phi_xi(N+1)           ! the basis and its derivative w.r.t. xi
   ! Local variables
   integer     :: i,j,m
   double precision :: tmp
   !
   ! Initialize variables
   phi      = 1.0d0
   phi_xi   = 0.0d0
   ! Lagrange polynomial and its derivative
   do m = 1, N+1
      do j = 1, N+1
         if (j .eq. m) cycle
         phi(m) = phi(m)*(xi-xGP(j))/(xGP(m)-xGP(j))
      end do

      do i = 1, N+1
         if( i .EQ. m) cycle
         tmp = 1.0d0
         do j = 1, N+1
            if( j .eq. i) cycle
            if( j .eq. m) cycle
            tmp = tmp*(xi-xGP(j))/(xGP(m)-xGP(j))
         end do
         phi_xi(m) = phi_xi(m) + tmp/(xGP(m)-xGP(i))
      end do
   end do

end subroutine BaseFunc1D
