function funcD_pe(Tsi)
    implicit NONE
    
     real*8 :: funcD_pe,Tsi,Tk
     real*8 :: Ep,rho_p,Rp,tau_g,a
     
     common/sil/Ep,rho_p,Rp,tau_g,a
     
     Tk = Tsi + 273.15d0
     funcD_pe = 22.86d0*a*Ep*Tk**0.5d0/tau_g
    
end function funcD_pe