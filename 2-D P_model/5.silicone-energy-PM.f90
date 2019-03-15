!Calculating "energy equation" of silicone
!Assume temperature is homogeneous.
subroutine Silicone_en_PM(j,k,Wsi,Wsi_after,Tsi,Tsi_after,Wair,Tair,Wair_after,Tair_after,Wsi_ave,Wrs_after,Tsi_after_new)
    implicit NONE
    
    integer :: Ni,Nj,Nk
    real*8 :: conv,del_z,del_rb,del_rp,del_t,start,finish
    real*8,dimension(Nj,Nk) :: Wair,Tair
    real*8,dimension(Nj,Nk) :: Wair_after,Tair_after,Wrs_after,Wsi_ave
    real*8,dimension(Ni,Nj,Nk) :: Wsi,Wsi_after,Tsi,Tsi_after,Tsi_after_new
    
    integer :: i,j,k
    real*8,dimension(Ni) :: C_w,C_p,C_e,RHS
    real*8 :: funcKd,funcCd,funcHads,funchc,funchm,funcCp_va
    real*8 :: Kd,Cd,Hads,hc,hm,Cp_va
    real*8 :: AA1,AA2,AA3,AA4,AA5,BB1,BB2,BB3,Bi,Ntu,DAR,gamma_b
    real*8 :: mu,rho,v,k_air,h_surr,Cp_v,T_surr
    real*8 :: Eb,rho_b,L,Rb,av
    real*8 :: Ep,rho_p,Rp,tau_g,a,tau,Re,Pi
    real*8 :: P_in,Tair_in,Wair_in
    
    common/int/Ni,Nj,Nk
    common conv,del_z,del_rb,del_rp,del_t
    common/sil/Ep,rho_p,Rp,tau_g,a
    common/flu/mu,rho,v,k_air,h_surr,Cp_v,T_surr
    common/bed/Eb,rho_b,L,Rb,av
    common/inlet/P_in,Tair_in,Wair_in
    common/other/tau,Re,Pi
    
    hc = funchc(Wair(j,k))
    hm = funchm(Tair(j,k))
    Hads = funcHads(Wsi(Ni,j,k))*1D3
    Cp_va = funcCp_va(Wair_after(j,k))
    Ntu = hm*P_in*L/v*4.d0/(2.d0*Rb)
    !Ntu = hm*P_in*L/v/(Pi*Rb**2)
    !Ntu = 22.65d0
    DAR = rho_b*L/rho/v/tau
    !DAR = 0.1285d0
    gamma_b = funcCd(Wsi_ave(j,k))/Cp_va
    
    !write(*,*)P_in,Ntu,DAR,gamma_b
    !pause
    
    BB1 = Ntu/DAR/gamma_b
    BB2 = hc/Cp_va/rho/hm
    BB3 = -Hads/Cp_va
    
    Tsi_after_new(:,j,k) = (BB1*BB2*Tair_after(j,k) + BB1*BB3*(Wrs_after(j,k) - Wair_after(j,k)) + Tsi(Ni,j,k)/del_t)/(1.d0/del_t + BB1*BB2)
    
    !do i = 1,Ni
    !    write(*,*)'k=',k,'j=',j,'i=',i,'Tsi',Tsi_after_new(i,j,k)
    !end do
    !pause
    !write(*,*)'k=',k,'j=',j,'Tsi',Tsi_after_new(1,j,k)
    !pause
    
end subroutine Silicone_en_PM