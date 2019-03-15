function funcWr_Wsi(Wsi,Tsi)
    implicit NONE
    
    real*8 :: Wsi,Tsi,funcWr_Wsi
    real*8 :: De,dRH,funcRH,funcPsat
    real*8 :: C0,C1,C2,C3,C4
    real*8 :: P_in,Tair_in,Wair_in,mu,rho,v,k_air,h_surr,Cp_v
    
    common/coeff/C0,C1,C2,C3,C4
    common/inlet/P_in,Tair_in,Wair_in
    common/flu/mu,rho,v,k_air,h_surr,Cp_v
    
    dRH = C1 + 2.d0*C2*Wsi + 3.d0*C3*Wsi**2.d0 + 4.d0*C4*Wsi**3.d0
    De = P_in - 0.378d0*funcRH(Wsi)*funcPsat(Tsi)
    
    funcWr_Wsi = (0.662d0*dRH*funcPsat(Tsi)*De - 0.378d0*dRH*funcPsat(Tsi)*0.662d0*funcRH(Wsi)*funcPsat(Tsi))/De**2
    !funcWr_Wsi = rho*dRH*funcPsat(Tsi)*0.662d0*P_in/De**2
    
end function funcWr_Wsi