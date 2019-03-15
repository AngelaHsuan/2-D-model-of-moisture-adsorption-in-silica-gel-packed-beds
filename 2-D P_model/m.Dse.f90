function funcD_se(Wsi,Tsi)
    implicit NONE
    
     real*8 :: funcD_se,Wsi,Tsi,Tk,funcHads
     
     Tk = Tsi + 273.15d0
     funcD_se = 1.6D-6*exp(-9.742D-4*funcHads(Wsi)*1D3/Tk)
    
end function funcD_se