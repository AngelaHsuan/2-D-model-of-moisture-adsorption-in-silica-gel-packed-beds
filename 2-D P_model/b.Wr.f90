function funcWr(Wsi,Tsi)
    implicit NONE
    
    real*8 :: funcWr,Wsi,Tsi,funcRH,funcPsat
    real*8 :: P_in
    common/inlet/P_in
    
    funcWr = (0.622d0*funcRH(Wsi) * funcPsat(Tsi))/(P_in - 0.378d0*funcRH(Wsi)*funcPsat(Tsi))
    
end function funcWr