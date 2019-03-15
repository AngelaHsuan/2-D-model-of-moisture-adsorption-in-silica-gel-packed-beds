function funchc(Wair)
    implicit NONE
    
     real*8 :: funchc,funcCp_va,Pr,Dp,Wair
     real*8 :: Ep,rho_p,Rp,mu,rho,v,k_air,tau,Re
     
     common/sil/Ep,rho_p,Rp
     common/flu/mu,rho,v,k_air
     common/other/tau,Re
     
     !Dp = 2.d0*Rp
     !Pr = mu*funcCp_va(Wair)/k_air
     !funchc = k_air/Dp*(2 + 1.1d0*Re**0.6d0*Pr**(1.d0/3.d0))
     funchc = 1.6d0*rho*v*Re**(-0.42d0)*funcCp_va(Wair)
    
end function funchc