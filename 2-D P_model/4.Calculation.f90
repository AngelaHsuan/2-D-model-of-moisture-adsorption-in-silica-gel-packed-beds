!
subroutine Calculation(Wsi_0,Tsi_0,Wrs_0,Wair_0,Tair_0,Wsi,Tsi,Wair,Tair)
    use omp_lib
    implicit NONE
    
    integer :: i,j,k,nn,kk,a,b,c
    integer :: Ni,Nj,Nk,Nk_t,dense,factor
    integer :: model
    real*8 :: err_k,err_j,tt,finish,start
    real*8 :: del_z,del_rp,del_rb,del_t,conv,w,del_zt
    real*8 :: Wsi_0,Tsi_0,Wrs_0,Wair_0,Tair_0
    real*8 :: P_in,Tair_in,Wair_in,W_all,W_all_0,W_tot
    real*8 :: Ep,rho_p,Rp,Eb,rho_b,L,Rb,tau,Re,Pi
    real*8 :: funcWr
    real*8,dimension(Nk) :: W_face_after,W_face
    real*8,dimension(Nj,Nk) :: Wsi_ave,Wsi_ave_after
    real*8,dimension(Nj,Nk) :: Wrs,Wrs_after,Wrs_after_new
    real*8,dimension(Nj,Nk) :: Wair,Tair,Wair_after,Tair_after,Wair_after_new,Tair_after_new
    real*8,dimension(Ni,Nj,Nk) :: Wsi,Tsi,Wsi_after,Tsi_after,Wsi_after_new,Tsi_after_new
    
    common/select/model
    common/int/Ni,Nj,Nk,Nk_t,dense,factor
    common conv,del_z,del_rb,del_rp,del_t,w,del_zt
    common/sil/Ep,rho_p,Rp
    common/bed/Eb,rho_b,L,Rb
    common/inlet/P_in,Tair_in,Wair_in
    common/other/tau,Re,Pi
    
    write(*,*)"Calculation starts!"
    !call cpu_time(start)
    !!$omp do private (i,j)
    !do k = 1,Nk
    !    do j = 1,Nj
    !        do i = 1,Ni
    !            Wsi(i,j,k) = Wsi_0
    !            Wsi_after(i,j,k) = Wsi_0
    !            Tsi(i,j,k) = Tsi_0
    !            Tsi_after(i,j,k) = Tsi_0
    !        end do
    !        Tair(j,k) = Tsi(Ni,j,k)
    !        Wsi_ave(j,k) = Wsi_0
    !        !Wface = Wsi_0
    !        Wrs(j,k) = Wrs_0
    !        Wair(j,k) = Wrs_0
    !
    !        Tair_after(j,k) = Tair(j,k)
    !        Wsi_ave_after(j,k) = Wsi_ave(j,k)
    !        Wrs_after(j,k) = Wrs(j,k)
    !        Wair_after(j,k) = Wair(j,k)
    !        end do
    !    end do
    !!$omp end do
    
    !Initial condition
    Wsi = Wsi_0
    Wsi_after = Wsi_0
    Tsi = Tsi_0
    Tsi_after = Tsi_0
    
    Tair = Tsi(Ni,:,:)
    Wsi_ave = Wsi_0
    !Wface = Wsi_0
    Wrs = Wrs_0
    Wair = Wrs_0
    
    Tair_after = Tair
    Wsi_ave_after = Wsi_ave
    Wrs_after = Wrs
    Wair_after = Wair
    !call cpu_time(finish)
    !write(*,*)finish-start
    
    W_all_0 = 0.d0
    do k = 1,Nk
        W_face(k) = Wair_0*(1.D0-Eb)*rho_p*del_z    !初始含水量
        W_all_0 = W_all_0 + W_face(k)
    end do
    
    !Boudary condition
    Tair(:,1)=Tair_in
    Wair(:,1)=Wair_in
    Tair_after(:,1)=Tair_in
    Wair_after(:,1)=Wair_in
    
    !Time proceeding
    do Nn = 1,1.d0/del_t,1
        !call cpu_time(start)
        tt = nn*del_t
        !write(*,'(A11,I8,A10,F15.6)')'timestep :',Nn,'time :',tt
        !write(*,*)''
        
        !z-step
        do kk = 1,Nk_t/dense*factor+1,1
            !write(*,*)nn,kk
            !call cpu_time(start)
            del_z = del_zt/factor
            err_k = 1D2
            do while (err_k > conv)
                err_k = 0.d0
                
                do j = 1,Nj
                    !write(*,*)'k',k
                    err_j = 0.001
                
                    !Doing iteration
                    do while (err_j > conv)   !Stop iteration while Wrs(j,k) converges
                        if (model == 1) then
                            call Silicone_ma_PM(j,kk,Wsi,Wsi_after,Tsi,Tsi_after,Wair,Wair_after,Tair,Tair_after,Wsi_ave,Wsi_ave_after,Wsi_after_new,Wrs_after)
                        else if (model == 2) then
                            call Silicone_ma(j,kk,Wsi,Wsi_after,Tsi,Tsi_after,Wair,Wair_after,Tair,Tair_after,Wsi_ave,Wsi_ave_after,Wsi_after_new,Wrs_after)
                        else
                            write(*,*)"mistake in 'calculation-1'"
                            stop
                        end if
                        Wsi_after(:,j,kk) = Wsi_after_new(:,j,kk)
                        
                        if (model == 1) then
                            call Silicone_en_PM(j,kk,Wsi,Wsi_after,Tsi,Tsi_after,Wair,Tair,Wair_after,Tair_after,Wsi_ave_after,Wrs_after,Tsi_after_new)
                        else if (model == 2) then
                            call Silicone_en(j,kk,Wsi,Wsi_after,Tsi,Tsi_after,Wair,Tair,Wair_after,Tair_after,Wsi_ave_after,Tsi_after_new)
                        else
                            write(*,*)"mistake in 'calculation-1'"
                            stop
                        end if
                        Tsi_after(:,j,kk) = Tsi_after_new(:,j,kk)
                        
                        !Calculate Wr and see that if it converges
                        Wrs_after_new(j,kk) = funcWr(Wsi_after(Ni,j,kk),Tsi_after(Ni,j,kk))
                        err_j = abs(Wrs_after_new(j,kk) - Wrs_after(j,kk))/abs(Wrs_after(j,kk))
                        Wrs_after(j,kk) = Wrs_after_new(j,kk)
                        !write(*,*)Wrs_after(j,kk)
                        !pause
                    end do  !while
                end do  !j
                
                call mass(kk,Wsi_ave_after,W_face_after)
                
                !call output(Wsi_after,Tsi_after,Wrs_after,Tair)
                !pause
            
                call Fluid_ma(kk,Wair,Tair,Wair_after,Tair_after,Wsi_after,Tsi_after,Wrs_after,Wair_after_new)
                call Fluid_en(kk,Wair,Tair,Wair_after,Tair_after,Wsi_after,Tsi_after,Wrs_after,Tair_after_new)
                
                !!$omp sections
                !!$omp section
                !call Fluid_ma(kk,Wair,Tair,Wair_after,Tair_after,Wsi_after,Tsi_after,Wrs_after,Wair_after_new)
                !!$omp section
                !call Fluid_en(kk,Wair,Tair,Wair_after,Tair_after,Wsi_after,Tsi_after,Wrs_after,Tair_after_new)
                !!$omp end sections
                
                do j = 1,Nj
                    if (err_k < (abs(Wair_after_new(j,kk) - Wair_after(j,kk))/abs(Wair_after(j,kk)))) then
                        err_k = abs(Wair_after_new(j,kk) - Wair_after(j,kk))/abs(Wair_after(j,kk))
                    end if
                    Wair_after(j,kk) = Wair_after_new(j,kk)
                    if (err_k < (abs(Tair_after_new(j,kk) - Tair_after(j,kk))/abs(Tair_after(j,kk)))) then
                        err_k = abs(Tair_after_new(j,kk) - Tair_after(j,kk))/abs(Tair_after(j,kk))
                    end if
                    Tair_after(j,kk) = Tair_after_new(j,kk)
                end do
                
            end do  !do while kk
            
            
            !call CPU_TIME(finish)
            !write(*,*)'time haha',finish - start
            !pause
        end do  !do kk
        
        do k = Nk_t/dense*factor+2,Nk,1
            del_z = del_zt
            err_k = 1D2
            do while (err_k > conv)
                err_k = 0.d0
                
                do j = 1,Nj
                    !write(*,*)'k',k
                    err_j = 0.001
                
                    !Doing iteration
                    do while (err_j > conv)   !Stop iteration while Wrs(j,k) converges
                        if (model == 1) then
                            call Silicone_ma_PM(j,k,Wsi,Wsi_after,Tsi,Tsi_after,Wair,Wair_after,Tair,Tair_after,Wsi_ave,Wsi_ave_after,Wsi_after_new,Wrs_after)
                        else if (model == 2) then
                            call Silicone_ma(j,k,Wsi,Wsi_after,Tsi,Tsi_after,Wair,Wair_after,Tair,Tair_after,Wsi_ave,Wsi_ave_after,Wsi_after_new,Wrs_after)
                        else
                            write(*,*)"mistake in 'calculation-2'"
                            stop
                        end if
                        Wsi_after(:,j,k) = Wsi_after_new(:,j,k)
                        
                        if (model == 1) then
                            call Silicone_en_PM(j,k,Wsi,Wsi_after,Tsi,Tsi_after,Wair,Tair,Wair_after,Tair_after,Wsi_ave_after,Wrs_after,Tsi_after_new)
                        else if (model == 2) then
                            call Silicone_en(j,k,Wsi,Wsi_after,Tsi,Tsi_after,Wair,Tair,Wair_after,Tair_after,Wsi_ave_after,Tsi_after_new)
                        else
                            write(*,*)"mistake in 'calculation-2'"
                            stop
                        end if
                        Tsi_after(:,j,k) = Tsi_after_new(:,j,k)
                        
                        !Calculate Wr and see that if it converges
                        Wrs_after_new(j,k) = funcWr(Wsi_after(Ni,j,k),Tsi_after(Ni,j,k))
                        err_j = abs(Wrs_after_new(j,k) - Wrs_after(j,k))/abs(Wrs_after(j,k))
                        Wrs_after(j,k) = Wrs_after_new(j,k)
                    
                    end do
                end do
                
                call mass(k,Wsi_ave_after,W_face_after)
                
                !call output(Wsi_after,Tsi_after,Wrs_after,Tair)
                !pause
                
                call Fluid_ma(k,Wair,Tair,Wair_after,Tair_after,Wsi_after,Tsi_after,Wrs_after,Wair_after_new)
                call Fluid_en(k,Wair,Tair,Wair_after,Tair_after,Wsi_after,Tsi_after,Wrs_after,Tair_after_new)
                !!$omp sections
                !!$omp section
                !call Fluid_ma(k,Wair,Tair,Wair_after,Tair_after,Wsi_after,Tsi_after,Wrs_after,Wair_after_new)
                !!$omp section
                !call Fluid_en(k,Wair,Tair,Wair_after,Tair_after,Wsi_after,Tsi_after,Wrs_after,Tair_after_new)
                !!$omp end sections
                
                !!$omp parallel sections
                !!$omp section
                !call Fluid_ma(k,Wair,Tair,Wair_after,Tair_after,Wsi_after,Tsi_after,Wrs_after,Wair_after_new)
                !!$omp section
                !call Fluid_en(k,Wair,Tair,Wair_after,Tair_after,Wsi_after,Tsi_after,Wrs_after,Tair_after_new)
                !!$omp end parallel sections
                
                do j = 1,Nj
                    if (err_k < (abs(Wair_after_new(j,k) - Wair_after(j,k))/abs(Wair_after(j,k)))) then
                        err_k = abs(Wair_after_new(j,k) - Wair_after(j,k))/abs(Wair_after(j,k))
                    end if
                    Wair_after(j,k) = Wair_after_new(j,k)
                    if (err_k < (abs(Tair_after_new(j,k) - Tair_after(j,k))/abs(Tair_after(j,k)))) then
                        err_k = abs(Tair_after_new(j,k) - Tair_after(j,k))/abs(Tair_after(j,k))
                    end if
                    Tair_after(j,k) = Tair_after_new(j,k)
                end do
                
            end do  !do while k
            
            !do j = 1,Nj
            !    write(*,*)'Wair k=',(k-Nk_t/dense*(factor-1)-1)*del_z,'j=',j,Wair_after(j,k)
            !end do
            !pause
            !do j = 1,Nj
            !    write(*,*)'Tair k=',(k-Nk_t/dense*(factor-1)-1)*del_z,'j=',j,Tair_after(j,k)
            !end do
            !pause
            
        end do  !do k
        
        
        !Refresh
        Wsi = Wsi_after
        Tsi = Tsi_after
        Wrs = Wrs_after
        Wsi_ave = Wsi_ave_after
        Wair = Wair_after
        Tair = Tair_after
        
        W_all = 0.d0
        do k = 1,Nk
            W_all = W_all + (W_face_after(k)*del_z)  !Wall:整個填充床的總吸水量
        end do
        W_tot = (W_all - W_all_0)*pi*Rb**2*L
        
        !if (mod(nn,50) == 0) then
        !    call animation(Wsi,Tsi,Wair,Tair,Wsi_ave)
        !end if
        
        if (mod(nn,60) == 0) then
            write(*,'(A11,I8,A10,F15.6)')'timestep :',Nn,'time :',tt
            !pause
        end if
        
        if (mod(nn,300) == 0) then
            !call output(Wsi,Tsi,Wair,Tair,Wsi_ave,tt)
            call output_sigmaplot(Tair,Tsi,Wair,Wsi,Wrs,Wsi_ave,W_tot,tt)
            !exit
            !call CPU_TIME(finish)
            !write(*,*)'time',finish - start
            !pause
        end if
        
        !if ((abs(Tair(1,Nk)-Tair_0)) <= 0.1d0 .AND. (abs(Wair(1,Nk)-Wair_in)) <= 1D-4 ) then    !!!!!!有改過收斂
        !    call output_sigmaplot(Tair,Tsi,Wair,Wsi,Wrs,Wsi_ave,W_tot,tt)
        !    write(*,*)'hi'
        !    exit
        !end if
        
        
    end do  !do nn
    
    !call output_sigmaplot(Tair,Tsi,Wair,Wsi,Wrs,Wsi_ave,W_tot,tt)
    
end subroutine Calculation