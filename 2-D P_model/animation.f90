!------------------------------------------------
!animation
subroutine animation(Wsi,Tsi,Wair,Tair,Wave)
    implicit NONE
    
    integer :: Ni,Nj,Nk,Nk_t,dense,factor
    integer :: i,j,k,kk
    real*8 :: conv,del_z,del_rb,del_rp,del_t
    real*8,dimension(Nj,Nk) :: Wair,Tair,Wave
    real*8,dimension(Ni,Nj,Nk) :: Wsi,Tsi
    
    integer :: Pad
    
    common/int/Ni,Nj,Nk,Nk_t,dense,factor
    common conv,del_z,del_rb,del_rp,del_t
    
    Pad = 100
    
    open (unit=Pad, file="animation-Tair.txt", status="UNKNOWN", ACCESS='APPEND')
    
    write(Pad,*) 'variables="z","rb","Tair"'
    write(Pad,*) 'zone i=', Nj,'j=',Nk,'DATAPACKING=POINT'
    
    del_z = 0.01d0/factor
    do kk = 1,Nk_t/dense*factor,1
        do j = 1,Nj
            write(Pad,'(3F20.12)') 0.d0+(kk-1)*del_z, 0.d0+(j-1)*del_rb, Tair(j,kk)
        end do
    end do
    del_z = 0.01d0
    do k = Nk_t/dense*factor+1,Nk,1
        do j = 1,Nj
            write(Pad,'(3F20.12)') 0.d0+(k-Nk_t/dense*(factor-1)-1)*del_z, 0.d0+(j-1)*del_rb, Tair(j,k)
        end do
    end do
    close (Pad,status = 'Keep')
    !---------------------------------------------------------------------------
    open (unit=Pad, file="animation-Wair.txt", status="UNKNOWN", ACCESS='APPEND')
    
    write(Pad,*) 'variables="z","rb","Wair"'
    write(Pad,*) 'zone i=', Nj,'j=',Nk,'DATAPACKING=POINT'
    
    del_z = 0.01d0/factor
    do kk = 1,Nk_t/dense*factor,1
        do j = 1,Nj
            write(Pad,'(3F20.12)') 0.d0+(kk-1)*del_z, 0.d0+(j-1)*del_rb, Wair(j,kk)
        end do
    end do
    del_z = 0.01d0
    do k = Nk_t/dense*factor+1,Nk,1
        do j = 1,Nj
            write(Pad,'(3F20.12)') 0.d0+(k-Nk_t/dense*(factor-1)-1)*del_z, 0.d0+(j-1)*del_rb, Wair(j,k)
        end do
    end do
    close (Pad,status = 'Keep')
    !---------------------------------------------------------------------------
    open (unit=Pad, file="animation-Wave.txt", status="UNKNOWN", ACCESS='APPEND')
    
    write(Pad,*) 'variables="z","rb","Wave"'
    write(Pad,*) 'zone i=', Nj,'j=',Nk,'DATAPACKING=POINT'
    
    del_z = 0.01d0/factor
    do kk = 1,Nk_t/dense*factor,1
        do j = 1,Nj
            write(Pad,'(3F20.12)') 0.d0+(kk-1)*del_z, 0.d0+(j-1)*del_rb, Wave(j,kk)
        end do
    end do
    del_z = 0.01d0
    do k = Nk_t/dense*factor+1,Nk,1
        do j = 1,Nj
            write(Pad,'(3F20.12)') 0.d0+(k-Nk_t/dense*(factor-1)-1)*del_z, 0.d0+(j-1)*del_rb, Wave(j,k)
        end do
    end do
    close (Pad,status = 'Keep')
    
end subroutine animation