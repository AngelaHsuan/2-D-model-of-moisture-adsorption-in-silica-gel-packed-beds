!Calculating "energy equation" of silicone
!Assume temperature distribution depends on r-dir only,
!so we do 1-D calculation.
subroutine Silicone_en(j,k,Wsi,Wsi_after,Tsi,Tsi_after,Wair,Tair,Wair_after,Tair_after,Wsi_ave,Tsi_after_new)
    implicit NONE
    
    integer :: Ni,Nj,Nk
    real*8 :: conv,del_z,del_rb,del_rp,del_t,start,finish
    real*8,dimension(Nj,Nk) :: Wsi_ave
    real*8,dimension(Nj,Nk) :: Wair,Tair,Wair_after,Tair_after,Wrs
    real*8,dimension(Ni,Nj,Nk) :: Wsi,Wsi_after,Tsi,Tsi_after,Tsi_after_new
    
    integer :: i,j,k
    real*8,dimension(Ni) :: C_w,C_p,C_e,RHS
    real*8 :: err,BC
    real*8 :: funcKd,funcCd,funcHads,funchc
    real*8 :: Kd,Cd,Hads,hc
    real*8 :: AA1,AA2,AA3,AA4,AA5,BB1,BB2,BB3,Bi
    real*8 :: Ep,rho_p,Rp,tau_g,a,tau,Re,Pi
    real*8 :: P_in,Tair_in,Wair_in
    character :: calculate
    character(len = 20) :: position
    
    common/int/Ni,Nj,Nk
    common conv,del_z,del_rb,del_rp,del_t
    common/sil/Ep,rho_p,Rp,tau_g,a
    common/inlet/P_in,Tair_in,Wair_in
    common/other/tau,Re,Pi
    
    calculate = 'E'
    
    hc = funchc(Wair(j,k))
    Cd = funcCd(Wsi_ave(j,k))
    
    do i = 1,Ni
        !Kd = funcKd(Wsi(i,j,k),Tsi(i,j,k))
        Hads = funcHads(Wsi(i,j,k))
        Kd = hc*Rp*10.d0
        BB1 = -rho_p*Hads*1D3/tau
        BB2 = rho_p*Cd/tau
        BB3 = Kd/Rp**2
        
        if (i == 1) then
        !Center of sphere
            position = 'center'
            call Silicone_coeff(i,j,k,AA1,AA2,AA3,AA4,AA5,BB1,BB2,BB3,Bi,BC,Wsi,Wsi_after,Tsi,Tsi_after,Wair,Tair,Wair_after,Tair_after,Wrs,calculate,position,C_w,C_p,C_e,RHS)
        else if (i == Ni) then
        !Surface of sphere
            !Bi = hc*Rp/Kd
            !write(*,*)Bi
            !pause
            Bi = 0.1d0
            position = 'surface'
            call Silicone_coeff(i,j,k,AA1,AA2,AA3,AA4,AA5,BB1,BB2,BB3,Bi,BC,Wsi,Wsi_after,Tsi,Tsi_after,Wair,Tair,Wair_after,Tair_after,Wrs,calculate,position,C_w,C_p,C_e,RHS)
        else
        !Other
            position = 'other'
            call Silicone_coeff(i,j,k,AA1,AA2,AA3,AA4,AA5,BB1,BB2,BB3,Bi,BC,Wsi,Wsi_after,Tsi,Tsi_after,Wair,Tair,Wair_after,Tair_after,Wrs,calculate,position,C_w,C_p,C_e,RHS)
        end if
    end do  !do i
    
    !write(*,*)'eeeeeeeeeeeeeeeeeeeeeen'
    !write(*,*)'k',k,'j',j
    !write(*,*)'C_w',C_w
    !write(*,*)'C_p',C_p
    !write(*,*)'C_e',C_e
    !write(*,*)'RHS',RHS
    !pause
    
    call TDMA(Ni,C_w,C_p,C_e,RHS,Tsi_after_new(:,j,k))
    
    !do i = 1,Ni
    !    write(*,*)'k=',k,'j=',j,'i=',i,'Tsi',Tsi_after_new(i,j,k)
    !end do
    !pause
    
end subroutine Silicone_en