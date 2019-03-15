function funchm(Tair)
    implicit NONE
    
     real*8 :: funchm,funcDair_H2O,Dp,Sc,Tair
     real*8 :: Ep,rho_p,Rp,mu,rho,v,k_air,tau,Re
     
     common/sil/Ep,rho_p,Rp
     common/flu/mu,rho,v,k_air
     common/other/tau,Re
     
     !Dp = 2.d0*Rp
     !Sc = mu/rho/funcDair_H2O(Tair)
     !funchm = funcDair_H2O(Tair)/Dp*(2.d0 + 1.1d0*Re**0.6d0*Sc**(1.d0/3.d0))
     funchm = 1.7d0*v*Re**(-0.42d0)
    
end function funchm