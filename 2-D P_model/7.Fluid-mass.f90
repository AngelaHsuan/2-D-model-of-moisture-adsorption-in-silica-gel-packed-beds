subroutine Fluid_ma(k,Wair,Tair,Wair_after,Tair_after,Wsi_after,Tsi_after,Wrs_after,Wair_after_new)
    implicit NONE
    
    integer :: Ni,Nj,Nk
    real*8 :: conv,del_z,del_rb,del_rp,del_t
    real*8,dimension(Nj,Nk) :: Wsi_ave,Wrs_after
    real*8,dimension(Nj,Nk) :: Wair,Tair,Wair_after,Tair_after,Wair_after_new
    real*8,dimension(Ni,Nj,Nk) :: Wsi_after,Tsi_after
    
    integer :: i,j,k
    real*8,dimension(Nj) :: C_s,C_p,C_n,RHS
    real*8 :: funchm,funcDair_H2O,err
    real*8 :: hm,Dair_H2O
    real*8 :: CC1,CC2,DD1,DD2
    real*8 :: Eb,rho_b,L,Rb,av,mu,rho,v,k_air,Cp_v,tau,Re,Pi
    real*8 :: P_in,Tair_in,Wair_in
    character :: calculate
    character(len = 20) :: position
    
    common/int/Ni,Nj,Nk
    common conv,del_z,del_rb,del_rp,del_t
    common/flu/mu,rho,v,k_air,Cp_v
    common/bed/Eb,rho_b,L,Rb,av
    common/other/tau,Re,Pi
    common/inlet/P_in,Tair_in,Wair_in
    
    calculate = 'M'
    
    if (k == 1) then    !Inlet
        Wair_after_new(:,k) = Wair_in
    else
        do j = 1,Nj
        !Calculate parameters
            CC1 = L*av*funchm(Tair_after(j,k))/Eb/v
            CC2 = L*funcDair_H2O(Tair_after(j,k))/Rb**2/Eb/v
            !CC2 = 0.d0 !for 1D
            if (j == Nj) then
            !Top of bed
                position = 'top'
                call Fluid_coeff(j,k,CC1,CC2,DD1,DD2,Wair_after,Tair_after,Wsi_after,Tsi_after,Wrs_after,calculate,position,C_s,C_p,C_n,RHS)
            else if (j == 1) then
            !Bottom of bed
                position = 'bottom'
                call Fluid_coeff(j,k,CC1,CC2,DD1,DD2,Wair_after,Tair_after,Wsi_after,Tsi_after,Wrs_after,calculate,position,C_s,C_p,C_n,RHS)
            else
            !Other
                position = 'other'
                call Fluid_coeff(j,k,CC1,CC2,DD1,DD2,Wair_after,Tair_after,Wsi_after,Tsi_after,Wrs_after,calculate,position,C_s,C_p,C_n,RHS)
            end if
        
        end do  !do j
        
        !if (k == Nk.OR.k == 50) then
        !write(*,*)'k',k,'j',j,'beeeeeeed'
        !write(*,*)'C_s',C_s
        !write(*,*)'C_p',C_p
        !write(*,*)'C_n',C_n
        !write(*,*)'RHS',RHS
        !pause
        !end if
        
        call TDMA(Nj,C_s,C_p,C_n,RHS,Wair_after_new(:,k))
            
        
    end if
    
    !err = maxval(abs(Wair_after_new-Wair_after))
    !Wair_after = Wair_after_new
    
    !if (k == Nk) then
    !    do j = 1,Nj
    !        write(*,*)'here k=',k,'j=',j,Wair_after(j,k)
    !    end do
    !    pause
    !end if
    
end subroutine Fluid_ma