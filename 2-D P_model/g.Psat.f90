function funcPsat(Tsi)
    implicit NONE
    
    real*8 :: funcPsat,Tsi,U,Tsi_F
    
    Tsi_F = Tsi*1.8d0 + 32.d0
    
    U = dexp(Tsi_F/100.D0)
    if (Tsi_F < 100.d0) then
        funcPsat = ((8.92222D-3*U + 0.247795d0)*U - 0.491179d0)*U + 0.274339d0
    else if (Tsi_F >= 100.d0 .and. Tsi_F < 150.d0) then
        funcPsat = ((-1.16464D-2*U + 0.411047d0)*U - 0.927337d0)*U + 0.66680d0
    else if (Tsi_F >= 150.d0 .and. Tsi_F < 200.d0) then
        funcPsat = ((-1.1018D-2*U + 0.396445d0)*U - 0.832936d0)*U + 0.48010d0
    else if (Tsi_F >= 200.d0) then
        funcPsat = ((-5.4575D-3*U + 0.269146d0)*U + 0.149028d0)*U - 2.0709d0
    else
        write(*,*)'Mistake in funcPsat!!'
        write(*,*)'Tsi_F=',Tsi_F
        pause
    end if
    
    
end function funcPsat