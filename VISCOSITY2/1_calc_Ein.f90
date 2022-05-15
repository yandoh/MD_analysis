!Copyright (c) 2022 Yoshimichi ANDOH
!Released under the MIT license
!https://opensource.org/licenses/mit-license.php

        use omp_lib
	implicit none
	integer(4)::i,j,k,il
        integer(4)::iam=0,nomp=1
        real(8),parameter :: dt=2.d-15
        real(8),parameter :: Temp=298.15d0     ! K
        integer(4),parameter :: iskip=100      ! t0のスキップ頻度
        integer(4),parameter::nfiles=50        ! 読み込みファイル数(配列allocateに利用)
        integer(4),parameter::nlperfile=500000 ! .mdmntrファイルの行数(MDステップ数)
	integer(4),parameter::nlines=nlperfile*nfiles
	real(8),parameter :: kB=1.380649d-23   ! Bolzman constant [J K−1]
!
	real(8) :: Pab(6,nlines)
        real(8) :: Pab0(6),Pab1(6)
	real(8) :: dG(6),dG2(6)
	real(8) :: avesPab(6,0:nlines)
        integer(4) :: nhit(0:nlines)
        real(8),allocatable :: avesPab_omp(:,:,:)
        integer(4),allocatable :: nhit_omp(:,:)
        integer(4)::nfilesin, nlinesin
	real(8) :: datin(20)
        real(8) :: kBT, volume, v_kBt

!$      nomp=omp_get_max_threads()
        allocate( avesPab_omp(6,0:nlines,0:nomp-1) )
        allocate( nhit_omp(0:nlines,0:nomp-1) )
        avesPab_omp=0d0
        nhit_omp=0

        kBT=kB*Temp ! J

	il=0
	volume=0d0
	avesPab=0d0
	nhit=0
        nfilesin=0

1	continue
        read(*,*,end=99) 
	read(*,*) 
	read(*,*) 
	read(*,*) 
	read(*,*) 
	read(*,*) 
        do j=1,nlperfile
          read(*,*,end=99) i, datin(1:20)
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
	v_kBt=volume/kBT   !  [m^3/J]
 	write(*,'(a,2es23.15,i10,i5,f10.3,i4)') "#",volume, v_kBt, &
     &   il, nfilesin, Temp, nomp

!! take difference
!$omp parallel default(none) &
!$omp private(iam,i,j,k,Pab1,dG,dG2) &
!$omp shared(il,Pab) &
!$omp shared(nhit_omp,avesPab_omp)
!$      iam=omp_get_thread_num()
!$omp do schedule(static,1)
 	do i=1,il-1,iskip
          dG=0d0
	  do j=i,il
            k=j-i
	    Pab1(1:6)=Pab(1:6,j)
            dG(1:6)=dG(1:6)+Pab1(1:6)*dt  !! tanzaku-sekibun
            dG2(1:6)=dG(1:6)**2
            nhit_omp(k,iam)=nhit_omp(k,iam)+1
            avesPab_omp(1:6,k,iam)=avesPab_omp(1:6,k,iam)+dG2(1:6)
	  enddo
	enddo
!$omp end do nowait
!$omp end parallel

!! summation over threads
        do k=0,nlines
        do iam=0,nomp-1
          nhit(k)=nhit(k)+nhit_omp(k,iam)
          avesPab(:,k)=avesPab(:,k)+avesPab_omp(:,k,iam)
        enddo
        enddo

!! take average & output
	do k=0,nlines-1
        if(k==0)then
          write(*,'(f18.3,6es23.15,i10)') k*dt*1e+12,  &
     &           0d0,0d0,0d0,0d0,0d0,0d0,nhit(k)
        else
          avesPab(1:6,k)=avesPab(1:6,k)/dble(nhit(k))
          write(*,'(f18.3,6es23.15,i10)') k*dt*1e+12,  &
     &           avesPab(1:6,k)*0.5d0*v_kBt, nhit(k)
        endif
	enddo

	stop
	end
