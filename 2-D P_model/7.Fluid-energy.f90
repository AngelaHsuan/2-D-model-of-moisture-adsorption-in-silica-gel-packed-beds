subroutine Fluid_en(k,Wair,Tair,Wair_after,Tair_after,Wsi_after,Tsi_after,Wrs_after,Tair_after_new)
    implicit NONE
    
    integer :: Ni,Nj,Nk
    real*8 :: conv,del_z,del_rb,del_rp,del_t
    real*8,dimension(Nj,Nk) :: Wsi_ave,Wrs_after
    real*8,dimension(Nj,Nk) :: Wair,Tair,Wair_after,Tair_after,Tair_after_new
    real*8,dimension(Ni,Nj,Nk) :: Wsi_after,Tsi_after
    
    integer :: i,j,k
    real*8,dimension(Nj) :: C_s,C_p,C_n,RHS
    real*8 :: mu,rho,v,k_air,Cp_v,Eb,rho_b,L,Rb,av,err
    real*8 :: P_in,Tair_in,Wair_in
    real*8 :: funchm,funcCp_va,funchc
    real*8 :: hm,Cp_va,hc
    real*8 :: CC1,CC2,DD1,DD2
    real*8 :: alpha,Kh,Km
    character :: calculate
    character(len = 20) :: position
    
    common/int/Ni,Nj,Nk
    common conv,del_z,del_rb,del_rp,del_t
    common/flu/mu,rho,v,k_air,Cp_v
    common/bed/Eb,rho_b,L,Rb,av
    common/inlet/P_in,Tair_in,Wair_in
    
    calculate = 'E'
    
    if (k == 1) then    !Inlet
        Tair_after_new(:,k) = Tair_in
    else
        do j = 1,Nj
            hm = funchm(Tair_after(j,k))
            Cp_va = funcCp_va(Wair_after(j,k))
            hc = funchc(Wair_after(j,k))
        !Calculate parameters
            alpha = k_air/rho/Cp_va
            Kh = L*av*hc/Eb/rho/v/Cp_va
            Km = L*av*hm/Eb/v
            DD1 = -Kh + Km*Cp_v/Cp_va*(Wair_after(j,k) - Wsi_after(Ni,j,k))
            DD2 = L*alpha/Rb**2/v
            !DD2 = 0.d0 !for 1D
            if (j == 1) then
            !Bottom of bed
                position = 'bottom'
                call Fluid_coeff(j,k,CC1,CC2,DD1,DD2,Wair_after,Tair_after,Wsi_after,Tsi_after,Wrs_after,calculate,position,C_s,C_p,C_n,RHS)
            else if (j == Nj) then
            !Top of bed
                position = 'top'
                call Fluid_coeff(j,k,CC1,CC2,DD1,DD2,Wair_after,Tair_after,Wsi_after,Tsi_after,Wrs_after,calculate,position,C_s,C_p,C_n,RHS)
            else
            !Other
                position = 'other'
                call Fluid_coeff(j,k,CC1,CC2,DD1,DD2,Wair_after,Tair_after,Wsi_after,Tsi_after,Wrs_after,calculate,position,C_s,C_p,C_n,RHS)
            end if
        end do  !do j
        
        !if (k == Nk.OR.k == 50) then
        !write(*,*)'k',k
        !write(*,*)'C_s',C_s
        !write(*,*)'C_p',C_p
        !write(*,*)'C_n',C_n
        !write(*,*)'RHS',RHS
        !pause
        !end if
        
        call TDMA(Nj,C_s,C_p,C_n,RHS,Tair_after_new(:,k))
            
    end if
        
        !do j = 1,Nj
        !    write(*,*)'Tair_after k=',k,'j=',j,Tair_after(j,k)
        !end do
        !pause
    
    !if (k == Nk) then
    !    do j = 1,Nj
    !        write(*,*)'k=',k,'j=',j,Tair_after(j,k)
    !    end do
    !    !pause
    !end if
    
end subroutine Fluid_en