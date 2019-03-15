!---------------------------------------------------------
! User
subroutine User
    implicit NONE
    
    integer :: model
    
    common/select/model
    
    write(*,*)'Which model do you want?'
    write(*,*)'(1) PM'
    write(*,*)'(2) K'
    read(*,*)model
    if (model .Ne. 1 .AND. model .Ne. 2) then
        write(*,*)"mistake in 'main'"
        stop
    end if
        
end subroutine User