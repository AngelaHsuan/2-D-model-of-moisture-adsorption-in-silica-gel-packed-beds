!Calculating "mass conservation" of silicone
!Assume temperature distribution depends on r-dir only,
!so we do 1-D calculation.
subroutine Silicone_ma(j,k,Wsi,Wsi_after,Tsi,Tsi_after,Wair,Wair_after,Tair,Tair_after,Wsi_ave,Wsi_ave_after,Wsi_after_new,Wrs)
    implicit NONE
    
    integer :: Ni,Nj,Nk
    real*8 :: conv,del_z,del_rb,del_rp,del_t,start,finish
    real*8,dimension(Nj,Nk) :: Wsi_ave,Wsi_ave_after
    real*8,dimension(Nj,Nk) :: Wair,Wair_after,Tair,Tair_after,Wrs
    real*8,dimension(Ni,Nj,Nk) :: Wsi,Wsi_after,Tsi,Tsi_after,Wsi_after_new
    
    integer :: i,j,k
    real*8,dimension(Ni) :: C_w,C_p,C_e,RHS
    real*8 :: err,sum,Wsi_i,rp_i
    real*8 :: Ep,rho_p,Rp,tau_g,a,mu,rho,v,k_air,Cp_v,tau,Re,Pi
    real*8 :: P_in,Tair_in,Wair_in
    real*8 :: funcKd,funcCd,funcHads,funcWr_Wsi,funchc,funcD_pe,funcWr_Tsi,funcD_se,funchm,funcPsat
    real*8 :: Kd,Cd,Hads,Wr_Wsi,hc,D_pe,Wr_Tsi,D_se,hm
    real*8 :: AA1,AA2,AA3,AA4,AA5,BB1,BB2,BB3,Bi,BC
    character :: calculate
    character(len = 20) :: position
    
    common/int/Ni,Nj,Nk
    common conv,del_z,del_rb,del_rp,del_t
    common/sil/Ep,rho_p,Rp,tau_g,a
    common/flu/mu,rho,v,k_air,Cp_v
    common/inlet/P_in,Tair_in,Wair_in
    common/other/tau,Re,Pi
    
    calculate = 'M'
    
    Cd = funcCd(Wsi_ave_after(j,k))
    hc = funchc(Wair(j,k))
    hm = funchm(Tair(j,k))
    
    do i = 1,Ni
        !Kd = funcKd(Wsi(i,j,k),Tsi(i,j,k))
        Kd = hc*Rp*10.d0
        Hads = funcHads(Wsi(i,j,k))*1D3
        Wr_Wsi = funcWr_Wsi(Wsi(i,j,k),Tsi(i,j,k))
        D_pe = funcD_pe(Tsi(i,j,k))
        Wr_Tsi = funcWr_Tsi(Wsi(i,j,k),Tsi(i,j,k))
        D_se = funcD_se(Wsi(i,j,k),Tsi(Ni,j,k))
        
        AA1 = rho/tau*(rho_p/rho + Ep*Wr_Wsi + D_pe*rho_p*Hads/Kd*Wr_Tsi)
        AA2 = Ep*rho/tau*(1.d0 - D_pe*rho_p*Cd/Ep/Kd)*Wr_Tsi
        AA3 = (D_se*rho_p + D_pe*rho*Wr_Wsi )/Rp**2
        
        if (i == 1) then
        !Center of sphere
            position = 'center'
            call Silicone_coeff(i,j,k,AA1,AA2,AA3,AA4,AA5,BB1,BB2,BB3,Bi,BC,Wsi,Wsi_after,Tsi,Tsi_after,Wair,Tair,Wair_after,Tair_after,Wrs,calculate,position,C_w,C_p,C_e,RHS)
        else if (i == Ni) then
        !Surface of sphere
            AA4 = -D_pe*rho/Rp**2/AA3*Wr_Tsi
            AA5 = rho*hm/Rp/AA3
            !Bi = hc*Rp/Kd
            Bi = 0.1d0
            !BC = rho*funchm(Tair(j,k))*Rp/(D_se+D_pe*Wr_Wsi/rho_p)/rho_p*(wair_after(j,k)-funcWr(Wsi(i,j,k),Tsi(i,j,k)))
            position = 'surface'
            call Silicone_coeff(i,j,k,AA1,AA2,AA3,AA4,AA5,BB1,BB2,BB3,Bi,BC,Wsi,Wsi_after,Tsi,Tsi_after,Wair,Tair,Wair_after,Tair_after,Wrs,calculate,position,C_w,C_p,C_e,RHS)
        else
        !Other
            position = 'other'
            call Silicone_coeff(i,j,k,AA1,AA2,AA3,AA4,AA5,BB1,BB2,BB3,Bi,BC,Wsi,Wsi_after,Tsi,Tsi_after,Wair,Tair,Wair_after,Tair_after,Wrs,calculate,position,C_w,C_p,C_e,RHS)
        end if                   
        !write(*,*)'ho',i,Tsi_after(i,j,k)
        !pause
        
    end do  !do i
    !!$omp end parallel do
    
    
    !write(*,*)'k',k,'j',j,'hiiiiiiiiiiiiiiiii'
    !write(*,*)'C_w',C_w
    !write(*,*)'C_p',C_p
    !write(*,*)'C_e',C_e
    !write(*,*)'RHS',RHS
    !pause
    
    !call TDMA(Ni,C_w,C_p,C_e,RHS,Wsi_after_i)
    !Wsi_after_new(:,j,k) = Wsi_after_i
    
    call TDMA(Ni,C_w,C_p,C_e,RHS,Wsi_after_new(:,j,k))
    
    !do i = 1,Ni
    !    write(*,*)'k=',k,'j=',j,'i=',i,'Wsi',Wsi_after_new(i,j,k)
    !end do
    !pause
    
    !Calculate Wsi_ave
    sum = 0.d0
    do i = 1,Ni
        rp_i = (i-1)*del_rp
        if (i == 1) then  !center
            Wsi_i = Wsi_after_new(i,j,k)*4.d0/3.d0*Pi*(del_rp/2.d0)**3.d0
        else if (i == Ni)then  !Surface
            Wsi_i = Wsi_after_new(i,j,k)*4.d0/3.d0*Pi*(1.5d0*(rp_i**2.d0)*del_rp - 0.75d0*rp_i*(del_rp**2.d0) + 0.125d0*(del_rp**3.d0))
        else
            Wsi_i = Wsi_after_new(i,j,k)*4.d0/3.d0*Pi*(3.d0*del_rp*(rp_i**2.d0) + 0.25d0*(del_rp**3.d0))
        end if
        sum = Wsi_i + sum
        
    end do
    Wsi_ave_after(j,k) = sum/((4.d0/3.d0)*Pi*(1.d0**3.d0))  !Radius is normalized
    
end subroutine Silicone_ma