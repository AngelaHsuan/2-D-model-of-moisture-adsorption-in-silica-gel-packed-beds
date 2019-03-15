!½Æ»s
program Main
    implicit NONE
    
    integer :: Ni,Nj,Nk
    integer :: model
    real*8 :: Wsi_0,Tsi_0,Wrs_0,Wair_0,Tair_0,P_in,Tair_in,Wair_in
    real*8,allocatable,dimension(:,:) :: Wair,Tair
    real*8,allocatable,dimension(:,:,:) :: Wsi,Tsi
    real*8 :: start,finish
    
    common/select/model
    common/int/Ni,Nj,Nk
    
    call cpu_time(start)
    call User
    
    if (model == 1) then
        !call Initial_K(Wsi_0,Tsi_0,Wrs_0,Wair_0,Tair_0)
        call Initial_PM(Wsi_0,Tsi_0,Wrs_0,Wair_0,Tair_0)
    else if (model == 2) then
        call Initial_K(Wsi_0,Tsi_0,Wrs_0,Wair_0,Tair_0)
        !call Initial_PM(Wsi_0,Tsi_0,Wrs_0,Wair_0,Tair_0)
    else
        write(*,*)"mistake in 'main'"
        stop
    end if
    
    allocate (Wsi(Ni,Nj,Nk),Tsi(Ni,Nj,Nk),Wair(Nj,Nk),Tair(Nj,Nk))
    
    call Calculation(Wsi_0,Tsi_0,Wrs_0,Wair_0,Tair_0,Wsi,Tsi,Wair,Tair)
    
    call cpu_time(finish)
    
    write(*,*)'time',finish - start
    
    
end program Main