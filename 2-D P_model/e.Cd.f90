function funcCd(Wsi_ave)
    implicit NONE
    
    real*8 :: Wsi_ave,funcCd
    
    funcCd = 4186.d0*Wsi_ave + 921.d0
    
end function funcCd