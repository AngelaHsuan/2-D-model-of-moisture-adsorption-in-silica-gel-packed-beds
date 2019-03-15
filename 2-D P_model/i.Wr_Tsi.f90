function funcWr_Tsi(Wsi,Tsi)
    implicit NONE
    
    real*8 :: funcWr_Tsi,funcWr,Wsi,Tsi,dT
     
    dT = 1D-6
    
    funcWr_Tsi = (funcWr(Wsi,Tsi+dT) - funcWr(Wsi,Tsi-dT))/2.d0/dT
    
end function funcWr_Tsi