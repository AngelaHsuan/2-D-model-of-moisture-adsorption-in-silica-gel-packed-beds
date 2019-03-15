function funcKd(Wsi,Tsi)
    implicit NONE
    
    real*8 :: Wsi,Tsi,funcKd
    
    funcKd = 0.37d0 + 0.97d0*Wsi + 0.0014d0*Tsi
    
end function funcKd