!---------------------------------------------------------
! Initial
subroutine Initial_PM(Wsi_0,Tsi_0,Wrs_0,Wair_0,Tair_0)
    implicit NONE
    
    integer :: Ni,Nj,Nk,Nk_t,dense,factor
    real*8 :: del_z,del_rp,del_rb,del_t,conv,w,del_zt
    real*8 :: Wrs
    real*8 :: Ep,rho_p,Rp,tau_g,a,tau_s,mu,rho
    real*8 :: P_in,Tair_in,Wair_in,v,k_air,h_surr,Cp_v,T_surr
    real*8 :: Eb,rho_b,L,Rb,av
    real*8 :: tau,Re,Pi
    real*8 :: C0,C1,C2,C3,C4
    real*8 :: Wsi_0,Tsi_0,Wair_0,Tair_0,Wrs_0
    real*8 :: funcWr
    real*8,allocatable,dimension(:) :: Wsi,Tsi
    real*8,allocatable,dimension(:,:) :: Wair,Tair
    
    common/int/Ni,Nj,Nk,Nk_t,dense,factor
    common conv,del_z,del_rb,del_rp,del_t,w,del_zt
    common/sil/Ep,rho_p,Rp,tau_g,a,tau_s
    common/flu/mu,rho,v,k_air,h_surr,Cp_v,T_surr
    common/bed/Eb,rho_b,L,Rb,av
    common/inlet/P_in,Tair_in,Wair_in
    common/other/tau,Re,Pi
    common/coeff/C0,C1,C2,C3,C4
    
    conv = 1D-4
    w = 1.d0
    dense = 2  !1/dense long of the input part having denser mesh
    factor = 1  !(factor) times denser
    del_zt = 0.01d0!0.01
    del_z = del_zt
    del_rp = 0.01d0!0.01
    !del_rp = 0.033d0
    del_rb = 0.01d0
    !del_t = 1.d0/240000.d0
    del_t = 1.d0/40000.d0
    !del_t = 1.d0/15000.d0
    
    Ni = 1.d0/del_rp + 1
    Nj = 1.d0/del_rb + 1
    Nk_t = 1.d0/del_z + 1
    Nk = Nk_t/dense*factor + Nk_t/dense*(dense-1) + 1
    !write(*,*)0.01d0/factor
    !pause
    
    !----------silicone properties------------
    
    Ep = 0.516d0    !溅ず场ふ回瞯
    rho_p = 1129.d0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Rp = 1.94D-3    !キА溅聋采畖
    tau_g = 2.8d0   !ふ回丁耎床Ρч
    tau_s = 2.8d0   !耎床Ρч
    a = 1.1D-9      !キА溅ふ回畖
    
    !------------fluid properties---------------
    
    mu = 1.9D-5
    rho = 1.2d0
    P_in = 14.696d0
    Tair_in = 23.3d0
    Wair_in = 0.01d0
    v = 0.21d0
    k_air = 0.0263d0
    !h_surr = 50.d0
    !k_air = 0.5d0
    !h_surr = 15000.d0
    h_surr = 0.d0
    T_surr = 20.d0
    !Dair_H2O = 2.79D-5
    Cp_v = 1884.d0
    
    !-------------bed properties--------------
    
    Eb = 0.6363d0
    rho_b = 721.1d0
    L = 0.0775d0
    Rb = 0.0388d0
    av = 3.d0*(1.d0-Eb)/Rp
    
    
    !-----Parameters of relative humidity-------
    
    C0 = 0.0078d0
    C1 = -0.05759d0
    C2 = 24.16554d0
    C3 = -124.4779d0
    C4 = 204.2264d0
    
    !--------Initial condition of silica gel------------
    
    Wsi_0 = 0.0417d0
    Tsi_0 = 23.3d0
    Wrs_0 = funcWr(Wsi_0,Tsi_0)
    
    !--------Initial condition of bed------------------
    
    Wair_0 = Wrs_0
    Tair_0 = 23.3d0
    
    !-----------Other system parameters----------
    
    tau = 25000.d0
    Re = 2.d0*rho*v*Rp/mu
    !Re = 49.30d0
    Pi = 3.14159265359d0
    
end subroutine Initial_PM