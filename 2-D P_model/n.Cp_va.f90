function funcCp_va(Wair)
    implicit NONE
    
     real*8 :: funcCp_va,Wair
     
     funcCp_va = 1884.d0*Wair + 1004.d0*(1.d0 - Wair)
    
end function funcCp_va