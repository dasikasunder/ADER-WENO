! boundary_conditions.f90
! author: sunder

!-----------------------------------------------------------------------
! Apply boundary conditions by filling ghost cell values
!-----------------------------------------------------------------------

subroutine apply_boundary_conditions
    use ader_weno
    implicit none
    
    integer :: ixx, oned_begin, oned_end, ilhs, irhs, iVar
    
    oned_begin = 1
    oned_end = IMAX
    
     ! Left end 
    
    do ixx = 1, nGhostCells
    
        ilhs = oned_begin - ixx
        
        if (bL .eq. 1) then ! Transmissive/Outflow/Continuitive boundary conditions
            irhs = oned_begin
        end if
        
        if (bL .eq. 3) then ! Periodic boundary conditions
            irhs = oned_end + 1 - ixx
        end if
        
        do iVar = 1, nVar
            uh(iVar,ilhs) = uh(iVar,irhs)
        end do
    
    end do

    ! Right end
    
    do ixx = 1, nGhostCells
        
        ilhs = oned_end + ixx
        
        if (bR .eq. 1) then ! Transmissive/Outflow/Continuitive boundary conditions
            irhs = oned_end
        end if
        
        if (bR .eq. 3) then ! Periodic boundary conditions
            irhs = oned_begin + ixx - 1
        end if
        
        do iVar = 1, nVar
            uh(iVar,ilhs) = uh(iVar,irhs)
        end do
    
    end do

end subroutine apply_boundary_conditions
