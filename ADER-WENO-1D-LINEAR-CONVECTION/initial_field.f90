! initial_field.f90
! author: sunder

!-----------------------------------------------------------------------
! Definition of initial condition
!-----------------------------------------------------------------------

subroutine initial_field(q0, xp)
    use ader_weno
    use constants, only : m_pi
    implicit none
    ! Argument list 
    double precision, intent(out) :: q0(nVar)
    double precision, intent(in)  :: xp
    ! Local variable declaration 
    double precision :: F1, G1, F2, G2, F3, G3
    double precision, parameter :: alpha = 10.0, a = 0.5, z = -0.7, delta = 0.005, beta = log(2.0)/(36.0*delta*delta)

    
    if (ICType .eq. 1) then
        q0(1) = sin(m_pi*xp)
    else if (ICType .eq. 2) then
        q0(1) = sin(m_pi*xp + (1.0d0/m_pi)*SIN(m_pi*xp))
    else if (ICType .eq. 3) then
        q0(1) = sin(m_pi*xp)**4
    else if (ICType .eq. 4) then
        if ((-0.5d0 .lt. xp) .AND. (xp .LT. 0.5d0)) then
            q0(1) = 1.0d0
        else
            q0(1) = 0.0d0
        end if
    else if (ICType .eq. 5) then
        F1 = sqrt(max(1.0 - alpha*alpha*(xp-a)*(xp-a) , 0.0d0 ))
        G1 = exp(-beta*(xp-z)*(xp-z))
        F2 = sqrt(max(1.0 - alpha*alpha*(xp-(a-delta))*(xp-(a-delta)) , 0.0d0 ))
        G2 = exp(-beta*(xp-(z-delta))*(xp-(z-delta)))
        F3 = sqrt(max(1.0 - alpha*alpha*(xp-(a+delta))*(xp-(a+delta)) , 0.0d0 ))
        G3 = exp(-beta*(xp-(z+delta))*(xp-(z+delta)))

        if (xp .ge. -0.8d0 .and. xp .lt. -0.6d0) then
            q0(1) = (1.0d0/6.0d0)*(4.0d0*G1 + G2 + G3)
        else if (xp .ge. -0.4d0 .and. xp .lt. -0.2d0) then
            q0(1) = 1.0d0
        else if (xp .ge. 0.0d0 .and. xp .lt. 0.2d0) then
            q0(1) = 1.0d0 - abs(10.0d0*(xp - 0.1d0))
        else if (xp .ge. 0.4d0 .and. xp .lt. 0.6d0) then
            q0(1) = (1.0d0/6.0d0)*(4.0d0*F1 + F2 + F3)
        else
            q0(1) = 0.0d0
        end if

    
    end if
    
end subroutine initial_field

