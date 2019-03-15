!----------------------------
subroutine Silicone_coeff(i,j,k,AA1,AA2,AA3,AA4,AA5,BB1,BB2,BB3,Bi,BC,Wsi,Wsi_after,Tsi,Tsi_after,Wair,Tair,Wair_after,Tair_after,Wrs,calculate,position,C_w,C_p,C_e,RHS)
    implicit NONE
    
    integer :: i,j,k,Ni,Nj,Nk
    real*8 :: conv,del_z,del_rb,del_rp,del_t
    real*8,dimension(Ni) :: C_w,C_p,C_e,RHS
    real*8,dimension(Nj,Nk) :: Wair,Tair,Wair_after,Tair_after,Wrs
    real*8,dimension(Ni,Nj,Nk) :: Wsi,Wsi_after,Tsi,Tsi_after
    real*8 :: funcWr
    character :: calculate
    character(len = 20) :: position
    
    real*8 :: AA1,AA2,AA3,AA4,AA5,BB1,BB2,BB3,Bi
    real*8 :: rp_i,rp_w,rp_e,BC
    
    common/int/Ni,Nj,Nk
    common conv,del_z,del_rb,del_rp,del_t
    
    rp_i = (i-1)*del_rp
    rp_w = (i-1.5d0)*del_rp  !rp-del_rp
    rp_e = (i-0.5d0)*del_rp  !rp+del_rp
    
    
    
    if (calculate == 'M') then
        if (position == 'center') then
            C_w(i) = 0.d0
            C_p(i) = AA1/del_t + 3.d0*8.d0*AA3*rp_e**2/(del_rp**4)
            C_e(i) = -3.d0*8.d0*AA3*rp_e**2/(del_rp**4)
            RHS(i) = AA1/del_t*Wsi(i,j,k) - AA2/del_t*(Tsi_after(i,j,k) - Tsi(i,j,k))
        else if (position == 'surface') then
            C_w(i) = -3.d0*AA3/del_rp*rp_w**2/(1.d0-rp_w**3)
            C_p(i) = AA1/del_t + 3.d0*AA3/del_rp*rp_w**2/(1.d0-rp_w**3)
            C_e(i) = 0.d0
            !BC = AA4*Bi*(Tair_after(j,k) - Tsi(Ni,j,k)) + AA5*(Wair_after(j,k) - Wrs(j,k))
            BC = AA4*Bi*(Tair_after(j,k) - Tsi_after(Ni,j,k)) + AA5*(Wair_after(j,k) - funcWr(Wsi_after(Ni,j,k),Tsi_after(Ni,j,k)))
            RHS(i) = AA1/del_t*Wsi(i,j,k) - AA2/del_t*(Tsi_after(i,j,k)-Tsi(i,j,k)) + 3.d0*AA3*BC*rp_i**2/(1.d0-rp_w**3.d0)
        else if (position == 'other') then
            C_w(i) = -3.d0*AA3*(rp_i/del_rp - 0.5d0)**2*del_rp/(rp_e**3-rp_w**3)
            C_p(i) = AA1/del_t + 3.d0*AA3*(2.d0*rp_i**2/del_rp**2 + 0.5d0)*del_rp/(rp_e**3-rp_w**3)
            C_e(i) = -3.d0*AA3*(rp_i/del_rp + 0.5d0)**2*del_rp/(rp_e**3-rp_w**3)
            RHS(i) = AA1/del_t*Wsi(i,j,k) - AA2/del_t*(Tsi_after(i,j,k) - Tsi(i,j,k))
            !write(*,*)AA1/del_t,3.d0*AA3*(2.d0*rp_i**2/del_rp**2 + 0.5d0)*del_rp,(rp_e**3-rp_w**3),del_rp
            !write(*,*)AA1/del_t, AA3*(2.d0*rp_i**2/del_rp**2 + 0.5d0)/(rp_i**2)
            !pause
            
            !AW = -3.D0*AA3*(RR2**2.D0)/RR3/dr
            !AP = AA1/dt + (3.D0*AA3*(RR1**2.D0+RR2**2.D0)/dr)/RR3
            !AE = -3.d0*AA3*(RR1**2.D0)/RR3/dr
            !SOURCE = AA1*W(I,N,K,1)/dt - AA2*(Trs(I,N,2)-Trs(I,N,1))/dt

            !write(*,*)'i',i
            !write(*,*)AA1/del_t
            !write(*,*)'C_w',C_w(i)
            !write(*,*)'C_p',C_p(i)
            !write(*,*)'C_e',C_e(i)
            !write(*,*)'RHS',RHS(i)
            !write(*,*)'hi',i,AA1,del_t,Wsi(i,j,k),AA2,Tsi_after(i,j,k),Tsi(i,j,k)
            !write(*,*)'ho',i,j,k,Tsi_after(i,j,k)
        else
            write(*,*)'Mistake : Silicone-coefficient -- M'
            stop
        end if
        
    else if (calculate == 'E') then
        if (position == 'center') then
            C_w(i) = 0.d0
            C_p(i) = BB2/del_t + 3.d0*8.d0*BB3*rp_e**2/(del_rp**4)
            C_e(i) = -3.d0*8.d0*BB3*rp_e**2/(del_rp**4)
            RHS(i) = BB2/del_t*Tsi(i,j,k) - BB1/del_t*(Wsi_after(i,j,k) - Wsi(i,j,k))
        else if (position == 'surface') then
            C_w(i) = -3.d0*BB3/del_rp*rp_w**2/(1.d0-rp_w**3)
            C_p(i) = BB2/del_t + 3.d0*BB3/del_rp*rp_w**2/(1.d0-rp_w**3)
            C_e(i) = 0.d0
            RHS(i) = BB2/del_t*Tsi(i,j,k) - BB1/del_t*(Wsi_after(i,j,k)-Wsi(i,j,k)) + 3.d0*BB3*(Bi*(Tair_after(j,k) - Tsi_after(Ni,j,k))*rp_i**2)/(1.d0-rp_w**3)
        else if (position == 'other') then
            C_w(i) = -3.d0*BB3*(rp_i/del_rp - 0.5d0)**2*del_rp/(rp_e**3-rp_w**3)
            C_p(i) = BB2/del_t + 3.d0*BB3*(2.d0*rp_i**2/del_rp**2 + 0.5d0)*del_rp/(rp_e**3-rp_w**3)
            C_e(i) = -3.d0*BB3*(rp_i/del_rp + 0.5d0)**2*del_rp/(rp_e**3-rp_w**3)
            RHS(i) = BB2/del_t*Tsi(i,j,k) - BB1/del_t*(Wsi_after(i,j,k) - Wsi(i,j,k))
        else
            write(*,*)'Mistake : Silicone-coefficient -- E'
            stop
        end if
    else
        write(*,*)'Mistake : Silicone-coefficient'
        stop
    end if
    
    
end subroutine Silicone_coeff