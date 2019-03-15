function funcDair_H2O(Tair)
    implicit NONE
    
    integer :: model
    real*8 :: funcDair_H2O,Tair,Tk
    real*8 :: P_in,Tair_in,Wair_in
    
    common/select/model
    common/inlet/P_in,Tair_in,Wair_in
    
    Tk = Tair + 273.15d0
    
    if (model == 1) then
        funcDair_H2O = 2.79D-5
    else if (model == 2) then
        funcDair_H2O = 1.735D-9*Tk**(1.685d0)/P_in
    else
        write(*,*)"mistake in 'funcDair_H2O'"
        stop
    end if
    
end function funcDair_H2O