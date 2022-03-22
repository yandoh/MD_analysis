!Copyright (c) 2022 Yoshimichi ANDOH
!Released under the MIT license
!https://opensource.org/licenses/mit-license.php

        use omp_lib
	implicit none
        integer(4)::nomp=1,iam=0
	integer(4)::i,j,k,il
        real(8),parameter :: Temp=298.15d0     ! 温度 [K] 
        integer(4),parameter :: iskip=100      ! t0のスキップ頻度
        integer(4),parameter::nfiles=50        ! 読み込みファイル数(配列allocateに利用)
        integer(4),parameter::nlperfile=500000 ! .mdmntrファイルの行数(MDステップ数)
	integer(4),parameter::nlines=nlperfile*nfiles
	real(8),parameter :: kB=1.380649d-23   ! Bolzman constant [J K−1] 
	real(8) :: Pab(6,nlines)
        real(8) :: Pab0(6),Pab1(6)
        real(8) :: Pab0xPab1(6)
        real(8) :: aveP01(6,0:nlines)
        real(8),allocatable :: aveP01_omp(:,:,:)
        integer(4) :: nhit(0:nlines)
        integer(4),allocatable :: nhit_omp(:,:)
        integer(4)::nfilesin, nlinesin
	real(8) :: datin(20)
        real(8) :: kBT, volume, v_kBt

!$      nomp=omp_get_max_threads()
        allocate( aveP01_omp(6,0:nlines,0:nomp-1) )
        allocate( nhit_omp(0:nlines,0:nomp-1) )
        aveP01_omp=0d0
        nhit=0

        kBT=kB*Temp ! J

	il=0
	volume=0d0
	nhit=0
	aveP01=0d0
        nfilesin=0

1	continue
        read(*,*,end=99) 
	read(*,*) 
	read(*,*) 
	read(*,*) 
	read(*,*) 
	read(*,*) 
        do j=1,nlperfile
	  read(*,*) i, datin(1:20)
	  il=il+1
	  volume=volume+datin(7)
	  Pab(1:6,il)=datin(15:20)
        enddo
        nfilesin=nfilesin+1
        if(nfilesin .gt. nfiles)then
          write(*,*) 'ERROR: nfilesin > nfiles, change nfiles value.'
          stop
        endif
	goto 1
99	continue

        volume=volume/dble(il)
	v_kBt=volume/kBT
 	write(*,'(a,2es23.15,i10,i5,f10.3,i4)') "#",volume, v_kBt, &
     &   il, nfilesin, Temp, nomp

        nlinesin=nfilesin*nlperfile

!! take autocorrelation
!$omp parallel default(none) &
!$omp private(iam,i,j,k,Pab0,Pab1,Pab0xPab1) &
!$omp shared(il,Pab) &
!$omp shared(nhit_omp,aveP01_omp)
!$      iam=omp_get_thread_num()
!$omp do schedule(static,1)
 	do i=1,il-1,iskip
	  Pab0(1:6)=Pab(1:6,i)
	  do j=i,il
	    Pab1(1:6)=Pab(1:6,j)
            Pab0xPab1(1:6)=Pab0(1:6)*Pab1(1:6)
            k=j-i
            nhit_omp(k,iam)=nhit_omp(k,iam)+1
            aveP01_omp(1:6,k,iam)=aveP01_omp(1:6,k,iam)+Pab0xPab1(1:6)
	  enddo
	enddo
!$omp end do nowait
!$omp end parallel

!! take sum over threads
	do k=0,nlinesin-1
        do iam=0,nomp-1
          nhit(k)=nhit(k)+nhit_omp(k,iam)
          aveP01(1:6,k)=aveP01(1:6,k)+aveP01_omp(1:6,k,iam)
        enddo
        enddo

!! take average
	do k=0,nlinesin-1
          aveP01(1:6,k)=aveP01(1:6,k)/dble(nhit(k))
          write(*,'(i10,6es23.15,i10)') k, aveP01(1:6,k)*v_kBt, nhit(k)
	enddo

	stop
	end
