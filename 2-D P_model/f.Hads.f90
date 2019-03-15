function funcHads(Wsi)
    implicit NONE
    
    integer :: model
    real*8 :: Wsi,funcHads
    
    common/select/model
    
    if (model == 1) then
        funcHads = -1079.D0*Wsi + 2745.D0
    else if (model == 2) then
        funcHads = -1079.D0*Wsi + 2745.D0
        !if (Wsi <= 0.05d0) then
        !    funcHads = -12400.d0*Wsi + 3500.d0
        !else if (Wsi > 0.05d0) then
        !    funcHads = -1400.d0*Wsi + 2950.d0
        !end if
    else
        write(*,*)"mistake in 'funcHads'"
        stop
    end if
    
end function funcHads