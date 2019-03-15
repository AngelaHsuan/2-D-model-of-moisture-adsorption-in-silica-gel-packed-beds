function funcRH(Wsi)
    implicit NONE
    
    real*8 :: C0,C1,C2,C3,C4
    real*8 :: Wsi,funcRH
    
    common/coeff/C0,C1,C2,C3,C4
    
    funcRH = C0 + C1*Wsi + C2*Wsi**2.d0 + C3*Wsi**3.d0 + C4*Wsi**4.d0
    
    
end function funcRH