subroutine mass(k,Wsi_ave_after,W_face_after)
    implicit NONE
    
    integer :: Ni,Nj,Nk,j,k
    real*8 :: W_face,sum
    real*8 :: A1,A2,A3,massin_t,massout_t,massin_t1,massout_t1,massin,massout,massin_tot,massout_tot,mass_tot
    real*8 :: conv,del_z,del_rb,mu,rho,v,Eb,rho_b,L,Rb,tau,Re,Pi,Ep,rho_p
    real*8,dimension(Nj,Nk) :: Wsi_ave_after
    real*8,dimension(Nk) :: W_face_after
    
    common/int/Ni,Nj,Nk
    common conv,del_z,del_rb
    common/sil/Ep,rho_p
    common/flu/mu,rho,v
    common/bed/Eb,rho_b,L,Rb
    common/other/tau,Re,Pi
    
    massout_t = 0.d0
    massin_t = 0.d0
    
    !do j = 1,Nj,1
    !    
    !    if (j == 1) then
    !        A1 = Pi*(0.5d0*del_rb)**2.d0 !係數在n=1
    !        massout_t1 = A1*Eb*Wair(j,Nk)*rho*v    !出口處每一個時間點*dt*period
    !        massin_t1 = A1*Eb*Wair(j,1)*rho*v          !進口處每一個時間點  
    !    else if (j == Nj) then
    !        A2 = Pi*(2.d0*(j-1)*del_rb*del_rb) !係數在n=2~NRMX
    !        massout_t1 = A3*Eb*Wair(j,Nk)*rho*v   !出口處每一個時間點
    !        massin_t1 = A3*Eb*Wair(j,1)*rho*v         !進口處每一個時間點
    !    else
    !        A3 = Pi*(((j-1)*del_rb*del_rb)-(0.25d0*(del_rb**2))) !係數在n=NRMXP1
    !        massout_t1 = A2*Eb*Wair(j,Nk)*rho*v   ! 出口處每一個時間點
    !        massin_t1 = A2*Eb*Wair(j,1)*rho*v         !進口處每一個時間點
    !    end if
    !    massout_t = massout_t1 + massout_t
    !    massin_t = massin_t1 + massin_t
    !end do
    !
    !massout = (massout_t/Pi)*Pi*(Rb**2.D0)       !算出平均後乘實際通過面積
    !massin = (massin_t/Pi)*Pi*(Rb**2.D0)
    !massin_tot = massin + massin_tot    !進口處累加水汽量(kg水汽)
    !massout_tot = massout + massout_tot !出口處累加水汽量(kg水汽) 
    
    
   !----
    sum = 0.d0
    do j = 1,Nj
        if (j == 1) then
            W_face = Wsi_ave_after(j,k)*pi*(0.5d0*del_rb)**2.d0
        else if (j == Nj) then
            W_face = Wsi_ave_after(j,k)*pi*(((j-1)*del_rb**2)-(0.25d0*(del_rb**2.d0)))
        else
            W_face = Wsi_ave_after(j,k)*pi*(2.D0*(j-1)*del_rb**2)
        end if
        sum = W_face + sum 
    end do
    W_face_after(k) = sum/pi*(1.d0-Eb)*rho_p    !為截面總吸水量
    
    
end