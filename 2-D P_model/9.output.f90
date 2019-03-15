!------------------------------------------------------------
! Result output
subroutine output(Wsi,Tsi,Wair,Tair,Wave,tt)
    implicit NONE
    
    integer :: Ni,Nj,Nk,Nk_t,dense,factor
    integer :: i,j,k,kk
    real*8 :: conv,del_z,del_rb,del_rp,del_t,tt
    real*8,dimension(Nj,Nk) :: Wair,Tair,Wave
    real*8,dimension(Ni,Nj,Nk) :: Wsi,Tsi
    
    integer :: Pad
    
    common/int/Ni,Nj,Nk,Nk_t,dense,factor
    common conv,del_z,del_rb,del_rp,del_t
    
    Pad = 100
    
    !Output Wsi
    !
    !open (unit=Pad, file="Wsi.txt", status="UNKNOWN")
    !
    !write(Pad,*) 'variables="rp","rb","z","Wsi"'
    !write(Pad,*) 'zone i=', Ni,'j=',Nj,'k=',Nk,'DATAPACKING=POINT'
    !
    !do k = 1,Nk,1
    !    do j = 1,Nj,1
    !        do i = 1,Ni,1
    !            write(Pad,'(4F20.12)') 0.d0+(i-1)*del_rp, 0.d0+(j-1)*del_rb, 0.d0+(k-1)*del_z, Wsi(i,j,k)
    !        end do
    !    end do
    !end do
    !
    !close (Pad,status = 'Keep')
    
    !!2D
    !Pad = 100
    !open (unit=Pad, file="Wsi.txt", status="UNKNOWN")
    !    
    !write(Pad,*) 'variables="rp","rb","Wsi"'
    !write(Pad,*) 'zone i=', Ni,'j=',Nj,'DATAPACKING=POINT'
    !
    !k = 1
    !do j = 1,Nj,1
    !    do i = 1,Ni,1
    !        write(Pad,'(4F20.12)') 0.d0+(i-1)*del_rp, 0.d0+(j-1)*del_rb, Wsi(i,j,k)
    !    end do
    !end do
    !
    !close (Pad,status = 'Keep')
    
    !Output Tsi
    !Pad = 100
    !
    !open (unit=Pad, file="Tsi.txt", status="UNKNOWN")
    !    
    !write(Pad,*) 'variables="rp","rb","z","Tsi"'
    !write(Pad,*) 'zone i=', Ni,'j=',Nj,'k=',Nk,'DATAPACKING=POINT'
    !
    !do k = 1,Nk,1
    !    do j = 1,Nj,1
    !        do i = 1,Ni,1
    !            write(Pad,'(4F20.12)') 0.d0+(i-1)*del_rp, 0.d0+(j-1)*del_rb, 0.d0+(k-1)*del_z, Tsi(i,j,k)
    !        end do
    !    end do
    !end do
    !
    !close (Pad,status = 'Keep')
    
    !2D
    !Pad = 100
    !
    !open (unit=Pad, file="Tsi.txt", status="UNKNOWN")
    !    
    !write(Pad,*) 'variables="rp","rb","Tsi"'
    !write(Pad,*) 'zone i=', Ni,'j=',Nj,'DATAPACKING=POINT'
    !
    !k = 1
    !do j = 1,Nj,1
    !    do i = 1,Ni,1
    !        write(Pad,'(4F20.12)') 0.d0+(i-1)*del_rp, 0.d0+(j-1)*del_rb, Tsi(i,j,k)
    !    end do
    !end do
    !
    !close (Pad,status = 'Keep')
    
    !Output Wair
    !Pad = 100
    !
    !open (unit=Pad, file="Wair.txt", status="UNKNOWN")
    !    
    !write(Pad,*) 'variables="z","rb","Wair"'
    !write(Pad,*) 'zone i=', Nk,'j=',Nj,'DATAPACKING=POINT'
    !
    !do j = 1,Nj,1
    !    do k = 1,Nk,1
    !        write(Pad,'(3F20.12)') 0.d0+(k-1)*del_z, 0.d0+(j-1)*del_rb, Wair(j,k)
    !    end do
    !end do
    !close (Pad,status = 'Keep')
    !
    !Output Tair
    !Pad = 100
    !
    !open (unit=Pad, file="Tair.txt", status="UNKNOWN")
    !    
    !write(Pad,*) 'variables="z","rb","Tair"'
    !write(Pad,*) 'zone i=', Nk,'j=',Nj,'DATAPACKING=POINT'
    !
    !do j = 1,Nj,1
    !    do k = 1,Nk,1
    !        write(Pad,'(3F20.12)') 0.d0+(k-1)*del_z, 0.d0+(j-1)*del_rb, Tair(j,k)
    !    end do
    !end do
    !close (Pad,status = 'Keep')
    
    !============================================================
    !1D
    Pad = 100
    
    open (unit=Pad, file="Wair1.txt", status="UNKNOWN")
        
    write(Pad,*) 'variables="z","Wair"'
    write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    
    del_z = 0.01d0/factor
    do kk = 1,Nk_t/dense*factor,1
        write(Pad,'(2F20.12)') 0.d0+(kk-1)*del_z, Wair(Nj,kk)
    end do
    del_z = 0.01d0
    do k = Nk_t/dense*factor+1,Nk,1
        write(Pad,'(2F20.12)') 0.d0+(k-Nk_t/dense*(factor-1)-1)*del_z, Wair(Nj,k)
    end do
    !write(Pad,*) 'variables="y","Wair"'
    !write(Pad,*) 'zone i=', Nj,'DATAPACKING=POINT'
    !do j = 1,Nj,1
    !    write(Pad,'(2F20.12)') 0.d0+(j-1)*del_rb, Wair(j,Nk)
    !end do
    
    close (Pad,status = 'Keep')
    
    
    !1D
    Pad = 100
    
    open (unit=Pad, file="Tair1.txt", status="UNKNOWN")
    
    write(Pad,*) 'variables="z","Tair"'
    write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    
    del_z = 0.01d0/factor
    do kk = 1,Nk_t/dense*factor,1
        write(Pad,'(2F20.12)') 0.d0+(kk-1)*del_z, Tair(Nj,kk)
    end do
    del_z = 0.01d0
    do k = Nk_t/dense*factor+1,Nk,1
        write(Pad,'(2F20.12)') 0.d0+(k-Nk_t/dense*(factor-1)-1)*del_z, Tair(Nj,k)
    end do
    !write(Pad,*) 'variables="y","Tair"'
    !write(Pad,*) 'zone i=', Nj,'DATAPACKING=POINT'
    !do j = 1,Nj,1
    !    write(Pad,'(2F20.12)') 0.d0+(j-1)*del_rb, Tair(j,Nk)
    !end do
    
    close (Pad,status = 'Keep')
    
    !Output Wave
    Pad = 100
    
    open (unit=Pad, file="Wave1.txt", status="UNKNOWN")
    
    write(Pad,*) 'variables="z","Wave"'
    write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    
    del_z = 0.01d0/factor
    do kk = 1,Nk_t/dense*factor,1
        write(Pad,'(2F20.12)') 0.d0+(kk-1)*del_z, Wave(Nj,kk)
    end do
    del_z = 0.01d0
    do k = Nk_t/dense*factor+1,Nk,1
        write(Pad,'(2F20.12)') 0.d0+(k-Nk_t/dense*(factor-1)-1)*del_z, Wave(Nj,k)
    end do
    !write(Pad,*) 'variables="y","Wave"'
    !write(Pad,*) 'zone i=', Nj,'DATAPACKING=POINT'
    !do j = 1,Nj,1
    !    write(Pad,'(2F20.12)') 0.d0+(j-1)*del_rb, Wave(j,Nk)
    !end do
    
    close (Pad,status = 'Keep')
    
    !==============================================================
    !Output Tair&Wair in 2D
    
    !open (unit=Pad, file="Tair-out.txt", status="UNKNOWN", ACCESS='APPEND')
    !write(Pad,'(2F20.12)') tt, Tair(1,Nk)
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Wair-out.txt", status="UNKNOWN", ACCESS='APPEND')
    !write(Pad,'(2F20.12)') tt, Wair(1,Nk)
    !close (Pad,status = 'Keep')
    
    
end subroutine output