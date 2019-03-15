!------------------------------------------------------------
! Result output for sigmaplot
subroutine output_sigmaplot(Tair,Tsi,Wair,Wsi,Wrs,Wave,W_tot,tt)
    implicit NONE
    
    integer :: Ni,Nj,Nk,Nk_t,dense,factor
    integer :: i,j,k,kk
    real*8 :: conv,del_z,del_rb,del_rp,del_t,tt,W_tot
    real*8,dimension(Nj,Nk) :: Wair,Tair,Wrs,Wave
    real*8,dimension(Ni,Nj,Nk) :: Wsi,Tsi
    
    integer :: Pad
    
    common/int/Ni,Nj,Nk,Nk_t,dense,factor
    common conv,del_z,del_rb,del_rp,del_t
    
    Pad = 100
    !Tair-----------------------------------------------
    !open (unit=Pad, file="Tair-j=1.txt", status="UNKNOWN")
    !    
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do k = 1,Nk,1
    !    write(Pad,'(1F20.12)')Tair(1,k)
    !end do
    !
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Tair-j=Nj.txt", status="UNKNOWN")
    !    
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do k = 1,Nk,1
    !    write(Pad,'(1F20.12)')Tair(Nj,k)
    !end do
    !
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Tair-j=Nj-half.txt", status="UNKNOWN")
    !    
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do k = 1,Nk,1
    !    write(Pad,'(1F20.12)')Tair(Nj/2,k)
    !end do
    !
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Tair-k=Nk-quarter.txt", status="UNKNOWN")
    !    
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do j = 1,Nj,1
    !    write(Pad,'(1F20.12)')Tair(j,Nk/4)
    !end do
    !
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Tair-k=Nk-half.txt", status="UNKNOWN")
    !    
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do j = 1,Nj,1
    !    write(Pad,'(1F20.12)')Tair(j,Nk/2)
    !end do
    !
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Tair-k=Nk-quarter3.txt", status="UNKNOWN")
    !    
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do j = 1,Nj,1
    !    write(Pad,'(1F20.12)')Tair(j,Nk*3/4)
    !end do
    !
    !close (Pad,status = 'Keep')
    !!Trs------------------------------------------------
    !!open (unit=Pad, file="Tsi.txt", status="UNKNOWN")
    !!    
    !!write(Pad,*) 'variables="z","Trs"'
    !!write(Pad,*) 'zone i=', Ni,'DATAPACKING=POINT'
    !!
    !!do i = 1,Ni,1
    !!    write(Pad,'(1F20.12)')Tsi(i,1,Nk)
    !!end do
    !!
    !!close (Pad,status = 'Keep')
    !!Wair------------------------------------------------
    !open (unit=Pad, file="Wair-j=1.txt", status="UNKNOWN")
    !
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do k = 1,Nk,1
    !    write(Pad,'(1F20.12)')Wair(1,k)
    !end do
    !
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Wair-j=Nj.txt", status="UNKNOWN")
    !    
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do k = 1,Nk,1
    !    write(Pad,'(1F20.12)')Wair(Nj,k)
    !end do
    !
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Wair-j=Nj-half.txt", status="UNKNOWN")
    !    
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do k = 1,Nk,1
    !    write(Pad,'(1F20.12)')Wair(Nj/2,k)
    !end do
    !
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Wair-k=Nk-quarter.txt", status="UNKNOWN")
    !    
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do j = 1,Nj,1
    !    write(Pad,'(1F20.12)')Wair(j,Nk/4)
    !end do
    !
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Wair-k=Nk-half.txt", status="UNKNOWN")
    !    
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do j = 1,Nj,1
    !    write(Pad,'(1F20.12)')Wair(j,Nk/2)
    !end do
    !
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Wair-k=Nk-quarter3.txt", status="UNKNOWN")
    !    
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do j = 1,Nj,1
    !    write(Pad,'(1F20.12)')Wair(j,Nk*3/4)
    !end do
    !
    !close (Pad,status = 'Keep')
    !!Wrs------------------------------------------------
    !open (unit=Pad, file="Wrs-j=1.txt", status="UNKNOWN")
    !
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do k = 1,Nk,1
    !    write(Pad,'(1F20.12)')Wrs(1,k)
    !end do
    !
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Wrs-j=Nj.txt", status="UNKNOWN")
    !    
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do k = 1,Nk,1
    !    write(Pad,'(1F20.12)')Wrs(Nj,k)
    !end do
    !
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Wrs-j=Nj-half.txt", status="UNKNOWN")
    !    
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do k = 1,Nk,1
    !    write(Pad,'(1F20.12)')Wrs(Nj/2,k)
    !end do
    !
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Wrs-k=Nk-quarter.txt", status="UNKNOWN")
    !    
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do j = 1,Nj,1
    !    write(Pad,'(1F20.12)')Wrs(j,Nk/4)
    !end do
    !
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Wrs-k=Nk-half.txt", status="UNKNOWN")
    !    
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do j = 1,Nj,1
    !    write(Pad,'(1F20.12)')Wrs(j,Nk/2)
    !end do
    !
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Wrs-k=Nk-quarter3.txt", status="UNKNOWN")
    !    
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do j = 1,Nj,1
    !    write(Pad,'(1F20.12)')Wrs(j,Nk*3/4)
    !end do
    !
    !close (Pad,status = 'Keep')
    !!==============================================================
    !!Wsi------------------------------------------------
    !open (unit=Pad, file="Wsi-j=1.txt", status="UNKNOWN")
    !
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do k = 1,Nk,1
    !    write(Pad,'(1F20.12)')Wsi(1,1,k)
    !end do
    !
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Wsi-j=Nj.txt", status="UNKNOWN")
    !    
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do k = 1,Nk,1
    !    write(Pad,'(1F20.12)')Wsi(1,Nj,k)
    !end do
    !
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Wsi-j=Nj-half.txt", status="UNKNOWN")
    !    
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do k = 1,Nk,1
    !    write(Pad,'(1F20.12)')Wsi(1,Nj/2,k)
    !end do
    !
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Wsi-k=Nk-quarter.txt", status="UNKNOWN")
    !    
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do j = 1,Nj,1
    !    write(Pad,'(1F20.12)')Wsi(1,j,Nk/4)
    !end do
    !
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Wsi-k=Nk-half.txt", status="UNKNOWN")
    !    
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do j = 1,Nj,1
    !    write(Pad,'(1F20.12)')Wsi(1,j,Nk/2)
    !end do
    !
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Wsi-k=Nk-quarter3.txt", status="UNKNOWN")
    !    
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do j = 1,Nj,1
    !    write(Pad,'(1F20.12)')Wsi(1,j,Nk*3/4)
    !end do
    !
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Wsi-i-1.txt", status="UNKNOWN")
    !    
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do i = 1,Ni,1
    !    write(Pad,'(1F20.12)')Wsi(i,Nj/2,Nk/4)
    !end do
    !
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Wsi-i-2.txt", status="UNKNOWN")
    !    
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do i = 1,Ni,1
    !    write(Pad,'(1F20.12)')Wsi(i,Nj/2,Nk/2)
    !end do
    !
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Wsi-i-3.txt", status="UNKNOWN")
    !    
    !!write(Pad,*) 'variables="z","Tair"'
    !!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do i = 1,Ni,1
    !    write(Pad,'(1F20.12)')Wsi(i,Nj/2,Nk*3/4)
    !end do
    !
    !close (Pad,status = 'Keep')
    
    !==============================================================
    !Output Tair&Wair in 2D
    
    !open (unit=Pad, file="Tair-out.txt", status="UNKNOWN", ACCESS='APPEND')
    !write(Pad,'(1F20.12)') Tair(1,Nk)
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Wair-out.txt", status="UNKNOWN", ACCESS='APPEND')
    !write(Pad,'(1F20.12)') Wair(1,Nk)
    !close (Pad,status = 'Keep')
    
    open (unit=Pad, file="mass-length.txt", status="UNKNOWN", ACCESS='APPEND')
    write(Pad,'(1F20.12)') W_tot
    close (Pad,status = 'Keep')
    
    !open (unit=Pad, file="Tair-1/4.txt", status="UNKNOWN", ACCESS='APPEND')
    !write(Pad,'(1F20.12)') Tair(1,Nk/4)
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Wair-1/4.txt", status="UNKNOWN", ACCESS='APPEND')
    !write(Pad,'(1F20.12)') Wair(1,Nk/4)
    !close (Pad,status = 'Keep')
    !
    !open (unit=Pad, file="Wave.txt", status="UNKNOWN", ACCESS='APPEND')
    !write(Pad,'(1F20.12)') Wave(1,Nk/2)
    !close (Pad,status = 'Keep')
    
    !==============================================================
    !1D
    !Pad = 100
    !
    !open (unit=Pad, file="Wair.txt", status="UNKNOWN")
    !    
    !write(Pad,*) 'variables="z","Wair"'
    !write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do kk = 1,Nk_t/dense*factor,1
    !    write(Pad,'(1F20.12)')Wair(1,kk)
    !end do
    !do k = Nk_t/dense*factor+1,Nk,1
    !    write(Pad,'(1F20.12)')Wair(1,k)
    !end do
    !!write(Pad,*) 'variables="y","Wair"'
    !!write(Pad,*) 'zone i=', Nj,'DATAPACKING=POINT'
    !!do j = 1,Nj,1
    !!    write(Pad,'(1F20.12)') Wair(j,Nk)
    !!end do
    !
    !close (Pad,status = 'Keep')
    !
    !
    !!1D
    !Pad = 100
    !
    !open (unit=Pad, file="Tair.txt", status="UNKNOWN")
    !
    !write(Pad,*) 'variables="z","Tair"'
    !write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !do kk = 1,Nk_t/dense*factor,1
    !    write(Pad,'(2F20.12)')Tair(1,kk)
    !end do
    !do k = Nk_t/dense*factor+1,Nk,1
    !    write(Pad,'(2F20.12)')Tair(1,k)
    !end do
    !!write(Pad,*) 'variables="y","Tair"'
    !!write(Pad,*) 'zone i=', Nj,'DATAPACKING=POINT'
    !!do j = 1,Nj,1
    !!    write(Pad,'(1F20.12)') Tair(j,Nk)
    !!end do
    !
    !close (Pad,status = 'Keep')
    !
    !!open (unit=Pad, file="Trs.txt", status="UNKNOWN")
    !!
    !!!write(Pad,*) 'variables="z","Trs"'
    !!!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !!!
    !!!do kk = 1,Nk_t/dense*factor,1
    !!!    write(Pad,'(2F20.12)')Tsi(Ni,Nj,kk)
    !!!end do
    !!!do k = Nk_t/dense*factor+1,Nk,1
    !!!    write(Pad,'(2F20.12)')Tsi(Ni,Nj,k)
    !!!end do
    !!write(Pad,*) 'variables="y","Tair"'
    !!write(Pad,*) 'zone i=', Nj,'DATAPACKING=POINT'
    !!do j = 1,Nj,1
    !!    write(Pad,'(2F20.12)') 0.d0+(j-1)*del_rb, Tair(j,Nk)
    !!end do
    !!
    !!close (Pad,status = 'Keep')
    !!
    !!open (unit=Pad, file="Wrs.txt", status="UNKNOWN")
    !!
    !!!write(Pad,*) 'variables="z","Wrs"'
    !!!write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !!!
    !!!do kk = 1,Nk_t/dense*factor,1
    !!!    write(Pad,'(2F20.12)')Wrs(Nj,kk)
    !!!end do
    !!!do k = Nk_t/dense*factor+1,Nk,1
    !!!    write(Pad,'(2F20.12)')Wrs(Nj,k)
    !!!end do
    !!write(Pad,*) 'variables="y","Tair"'
    !!write(Pad,*) 'zone i=', Nj,'DATAPACKING=POINT'
    !!do j = 1,Nj,1
    !!    write(Pad,'(2F20.12)') 0.d0+(j-1)*del_rb, Tair(j,Nk)
    !!end do
    !!
    !!close (Pad,status = 'Keep')
    !
    !
    !open (unit=Pad, file="Wave.txt", status="UNKNOWN")
    !
    !write(Pad,*) 'variables="z","Wave"'
    !write(Pad,*) 'zone i=', Nk,'DATAPACKING=POINT'
    !
    !del_z = 0.01d0/factor
    !do kk = 1,Nk_t/dense*factor,1
    !    write(Pad,'(F20.12)')  Wave(1,kk)
    !end do
    !del_z = 0.01d0
    !do k = Nk_t/dense*factor+1,Nk,1
    !    write(Pad,'(F20.12)')Wave(1,k)
    !end do
    !!write(Pad,*) 'variables="y","Wave"'
    !!write(Pad,*) 'zone i=', Nj,'DATAPACKING=POINT'
    !!do j = 1,Nj,1
    !!    write(Pad,'(1F20.12)') Wave(j,Nk)
    !!end do
    !
    !close (Pad,status = 'Keep')
    
    
    
    
    
end subroutine output_sigmaplot