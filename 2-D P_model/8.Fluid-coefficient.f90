!----------------------------
subroutine Fluid_coeff(j,k,CC1,CC2,DD1,DD2,Wair_after,Tair_after,Wsi_after,Tsi_after,Wrs_after,calculate,position,C_s,C_p,C_n,RHS)
    implicit NONE
    
    integer :: j,k,Ni,Nj,Nk
    real*8 :: conv,del_z,del_rb,del_rp,del_t
    real*8 :: mu,rho,v,k_air,h_surr,Cp_v,T_surr,Eb,rho_b,L,Rb,av,tau,Re,Pi
    real*8,dimension(Nj) :: C_s,C_p,C_n,RHS
    real*8,dimension(Nj,Nk) :: Wair_after,Tair_after,Wrs_after
    real*8,dimension(Ni,Nj,Nk) :: Wsi_after,Tsi_after
    character :: calculate
    character(len = 20) :: position
    
    real*8 :: CC1,CC2,DD1,DD2,BC
    real*8 :: rb_j,rb_s,rb_n
    
    common/int/Ni,Nj,Nk
    common conv,del_z,del_rb,del_rp,del_t
    common/flu/mu,rho,v,k_air,h_surr,Cp_v,T_surr
    common/bed/Eb,rho_b,L,Rb,av
    common/other/tau,Re,Pi
    
    rb_j = (j-1)*del_rb
    rb_n = (j-0.5d0)*del_rb
    rb_s = (j-1.5d0)*del_rb
    
    if (calculate == 'M') then
        if (position == 'bottom') then
            C_s(j) = 0.d0
            C_p(j) = 0.125d0*CC1*del_rb + CC2*rb_n/del_rb**2 + 0.25d0*del_rb/del_z
            C_n(j) = -CC2*rb_n/del_rb**2
            RHS(j) = (0.25d0*del_rb/del_z - 0.125d0*CC1*del_rb - CC2*rb_j/del_rb**2)*Wair_after(j,k-1)&
                    & + 0.125d0*CC1*del_rb*(Wrs_after(j,k) + Wrs_after(j,k-1)) + CC2*rb_j/del_rb**2*Wair_after(j+1,k-1)
        else if (position == 'top') then
            C_s(j) = -CC2*rb_s/del_rb**2
            C_p(j) = 0.5d0*CC1*(1.d0-0.25d0*del_rb) + CC2*rb_s/del_rb**2 + (1.d0-0.25d0*del_rb)/del_z
            C_n(j) = 0.d0
            RHS(j) = ((1.d0-0.25d0*del_rb)/del_z - 0.5d0*CC1*(1.d0-0.25d0*del_rb) - CC2*rb_s/del_rb**2)*Wair_after(j,k-1)&
                    & + 0.5d0*CC1*(1.d0-0.25d0*del_rb)*(Wrs_after(j,k) + Wrs_after(j,k-1)) + CC2*rb_s/del_rb**2*Wair_after(j-1,k-1)
        else if (position == 'other') then
            C_s(j) = -0.5d0*CC2*rb_s/rb_j/del_rb**2
            C_p(j) = 0.5d0*CC1 + CC2/del_rb**2 + 1.d0/del_z
            C_n(j) = -0.5d0*CC2*rb_n/rb_j/del_rb**2
            RHS(j) = (1.d0/del_z - 0.5d0*CC1 - CC2/del_rb**2)*Wair_after(j,k-1) + 0.5d0*CC1*(Wrs_after(j,k) + Wrs_after(j,k-1))&
                    & + 0.5d0*CC2*rb_n/del_rb**2/rb_j*Wair_after(j+1,k-1) + 0.5d0*CC2*rb_s/del_rb**2/rb_j*Wair_after(j-1,k-1)
        else
            write(*,*)'Mistake : Fluid-coefficient -- M'
            stop
        end if
    else if (calculate == 'E') then
        if (position == 'bottom') then
            C_s(j) = 0.d0
            C_p(j) = -0.125d0*DD1*del_rb + DD2*rb_n/del_rb**2 + 0.25d0*del_rb/del_z
            C_n(j) = -DD2*rb_n/del_rb**2
            RHS(j) = (0.25d0*del_rb/del_z + 0.125d0*DD1*del_rb - DD2*rb_j/del_rb**2)*Tair_after(j,k-1)&
                    & - 0.125d0*DD1*del_rb*(Tsi_after(Ni,j,k) + Tsi_after(Ni,j,k-1)) + DD2*rb_j/del_rb**2*Tair_after(j+1,k-1)
        else if (position == 'top') then
            !alpha = k_air/rho/funcCp_va(Wair_after(j,k))
            !Kh = L*av*funchc(Wair_after(j,k))/Eb/rho/v/funcCp_va(Wair_after(j,k))
            !Km = L*av*funchm(Tair_after(j,k))/Eb/v
            !DD1 = -Kh + Km*Cp_v/funcCp_va(Wair_after(j,k))*(Wair_after(j,k) - Wsi_after(Ni,j,k))
            !DD2 = L*alpha/Rb**2/v
            BC = (-Rb*h_surr/k_air)
            C_s(j) = -DD2*rb_s/del_rb**2
            C_p(j) = -0.5d0*DD1*(1.d0-0.25d0*del_rb) + DD2*rb_s/del_rb**2 + (1.d0-0.25d0*del_rb)/del_z - DD2*BC/del_rb
            C_n(j) = 0.d0
            RHS(j) = ((1.d0-0.25d0*del_rb)/del_z + 0.5d0*DD1*(1.d0-0.25d0*del_rb) - DD2*rb_s/del_rb**2 + DD2*BC/del_rb)*Tair_after(j,k-1)&
                    & - 0.5d0*DD1*(1.d0-0.25d0*del_rb)*(Tsi_after(Ni,j,k) + Tsi_after(Ni,j,k-1)) + DD2*rb_s/del_rb**2*Tair_after(j-1,k-1)&
                    & - 2.d0*DD2*BC/del_rb*T_surr
            !BC = (-Rb*h_surr/k_air)
            !C_s(j) = -2.d0*DD2*rb_s/del_rb**2
            !C_p(j) = -DD1*(1.d0-0.25d0*del_rb) + 2.d0*DD2*rb_s/del_rb**2 - 2.d0*DD2*BC/del_rb + (1.d0-0.25d0*del_rb)/del_z
            !C_n(j) = 0.d0
            !RHS(j) = (1.d0-0.25d0*del_rb)/del_z*Tair_after(j,k-1) - DD1*(1.d0-0.25d0*del_rb)*Tsi_after(Ni,j,k) - 2.d0*DD2*BC/del_rb*T_surr
            
            !TT1=2.D0*PI*DD2
            !TT2=2.D0*PI*DD2*hsurr/k_air
            !rn=0.d0
            !rp=(N-1.D0)*dn
            !rs=(N-2.D0)*dn
            !RR1=(rp+rn)/2.d0
            !RR2=(rs+rp)/2.D0
            !TAS=-2.D0*PI*DD2*rb_s/del_rb
            !TAP=1.d0/del_z+(2.D0*PI*DD2*rb_s/del_rb)-(2.D0*PI*DD2*BC/Rb/del_z) - DD1
            !TAN=0.D0
            !TSOURCE=Tair(I-1,N,2)/del_z-(2.D0*PI*DD2*BC/Rb*Tsurr/del_z)-DD1*Tsi(I,N,KRMXP1,2)
            
            
        else if (position == 'other') then
            C_s(j) = -0.5d0*DD2*rb_s/rb_j/del_rb**2
            C_p(j) = -0.5d0*DD1 + DD2/del_rb**2 + 1.d0/del_z
            C_n(j) = -0.5d0*DD2*rb_n/rb_j/del_rb**2
            RHS(j) = (1.d0/del_z + 0.5d0*DD1 - DD2/del_rb**2)*Tair_after(j,k-1) - 0.5d0*DD1*(Tsi_after(Ni,j,k) + Tsi_after(Ni,j,k-1))&
                    & + 0.5d0*DD2*rb_n/del_rb**2/rb_j*Tair_after(j+1,k-1) + 0.5d0*DD2*rb_s/del_rb**2/rb_j*Tair_after(j-1,k-1)
        else
            write(*,*)'Mistake : Fluid-coefficient -- E'
            stop
        end if
    else
        write(*,*)'Mistake : Fluid-coefficient'
        stop
    end if
    
    
end subroutine Fluid_coeff