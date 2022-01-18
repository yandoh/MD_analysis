!Copyright (c) 2022 Yoshimichi ANDOH
!Released under the MIT license
!https://opensource.org/licenses/mit-license.php

! === inputs ===
! (1) ./sysinfo
! (2) ./massinfo.mdff
! (3) ./dcd
!
        implicit none
        integer(4),parameter::ncopy=1
!
	integer(4),allocatable::nav(:),nmv(:)
	integer(4)::h,i,j,k,l,m,n,komtot,kom,ix,iy,iz
	real(8)::box(3)
	real(8)::alpha,beta,gamma
	real(8),allocatable :: mass(:,:), molmass(:)
	real(8)::avec(3),bvec(3),cvec(3)
	integer(4) :: totmol
	integer(4) :: totnp
	real(8) :: smass
      real(8) :: HH(3,3), rH(3,3), tHHH(3,3)
	real(8),allocatable :: pos(:,:,:)
	real(8),allocatable :: cntmol(:,:,:)
	real(8),allocatable :: dvecz(:,:)
      real(8) :: cntlyr(2,3)
	real(8) :: lx,ly,lz
	real(8) :: piQ
	real(8) :: cutr
!dcd
      character(len=4)::Aname
      integer(4) :: nstr, nptot, iflame, nflame
      integer(4),dimension(20)::icntrl
      real(8),dimension(6)::cellstr
      real(4),allocatable::flpx(:),flpy(:),flpz(:)
      real(8)::pi
      Aname='CORD'
      icntrl=0
      nstr=0
      pi=acos(-1d0)
!dcd

	pi=dacos(-1d0)
	piQ=pi/4d0

!### input sysinfo ###
	open(21,file='sys_info')
        read(21,*) ! skip
 	read(21,*) komtot
 	allocate( nav(komtot), nmv(komtot) )
        do kom=1,komtot
          read(21,*) nav(kom),nmv(kom)
        enddo
	close(21)

	totmol=0
	totnp=0
	do kom=1,komtot
	  totmol=totmol+nmv(kom)
	  totnp =totnp +nmv(kom)*nav(kom)
	enddo
	allocate( mass(komtot,totnp), molmass(komtot) )
	allocate( pos(3,komtot,totnp) )
	allocate( cntmol(3,komtot,totmol) )
	allocate( dvecz(komtot,totmol) )
	
!### input mass ###
	open(3,file='massinfo.mdff')
	do kom=1,komtot
	read(3,*) ! skip np
	do m=1,nav(kom)!*nmv(kom)
	  read(3,*) mass(kom,m)
	enddo
	read(3,*) ! skip np
	enddo
	close(3)

	molmass=0d0
	do kom=1,komtot
	do m=1,nav(kom)
	  molmass(kom)=molmass(kom)+mass(kom,m)
	enddo
	enddo

!### input .dcd file ###
	open(34,file='./DCD',form='unformatted')
	read(34) Aname,icntrl
	read(34) nstr
	read(34) nptot
      allocate( flpx(nptot),flpy(nptot),flpz(nptot) )

	nflame=icntrl(1)
	write(*,'(a,i7,2i5)') '#', nflame,totmol,(ncopy+1+ncopy)

	DO iflame=1,icntrl(1)
!debug
!!!	write(9,'(a,i10)') '#', iflame
!debug
!!	if(mod(iflame,100)==0) write(*,*) 'Converting ',iflame,'-th flame'
	read(34) cellstr !*1e+10
	read(34) (flpx(i),i=1,nptot)
	read(34) (flpy(i),i=1,nptot)
	read(34) (flpz(i),i=1,nptot)

	k=0
	do kom=1,komtot
	do h=1,nmv(kom)
	i=(h-1)*nav(kom)
	  do m=1,nav(kom)
	k=k+1
          pos(1,kom,i+m)=flpx(k)
          pos(2,kom,i+m)=flpy(k)
          pos(3,kom,i+m)=flpz(k)
        enddo
      enddo
	enddo
 
      box(1)=cellstr(1)   ! |a|
      gamma =cellstr(2)   ! degree
      box(2)=cellstr(3)   ! |b|
      beta  =cellstr(4)   ! degree
      alpha =cellstr(5)   ! degree
      box(3)=cellstr(6)   ! |c|

      alpha=alpha/180d0*pi
      beta =beta /180d0*pi
      gamma=gamma/180d0*pi

      call creatematrix(  &
                box,alpha,beta,gamma,HH,rH,tHHH)
	avec(:)=HH(:,1)
	bvec(:)=HH(:,2)
	cvec(:)=HH(:,3)

!### calc. COM ###
	cntmol=0d0
	DO KOM=1,komtot
	do h=1,nmv(KOM)
	i=(h-1)*nav(KOM)
	do k=1,3  ! x, y, z
	  do m=1,nav(KOM)
	    cntmol(k,kom,h)=cntmol(k,kom,h)+pos(k,kom,i+m)*mass(kom,m)
	  enddo
	  cntmol(k,kom,h)=cntmol(k,kom,h)/molmass(kom)
	enddo
	enddo
	ENDDO

!	^^^ output cntmol ^^^
 	do kom=1,komtot
 	do ix=-ncopy,+ncopy ; do iy=-ncopy,+ncopy ; do iz=-ncopy,+ncopy
 	do h=1,nmv(kom)
	  write(*,'(3f12.5,i3)') cntmol(:,kom,h) &
     &                           +ix*avec(:) + iy*bvec(:) + iz*cvec(:)
 	enddo ! h
	enddo; enddo; enddo
 	enddo ! kom

 	write(*,'(a1,3f12.5)') '#',alpha,beta,gamma
 	write(*,'(a1,3f12.5)') '#',box

	ENDDO
 	close(34)

	stop
	end

!########################################################################
        subroutine reverse(HH,rH)
!########################################################################
        real(8)::HH(3,3),rH(3,3)
        real(8)::det,rdet

        det=HH(1,1)*( HH(2,2)*HH(3,3)-HH(2,3)*HH(3,2) )  &
           -HH(1,2)*( HH(2,1)*HH(3,3)-HH(2,3)*HH(3,1) )  &
           +HH(1,3)*( HH(2,1)*HH(3,2)-HH(2,2)*HH(3,1) )  
        if(det.eq.0d0) stop 'det = 0'
        rdet=1d0/det

        rH(1,1)=rdet*( HH(2,2)*HH(3,3)-HH(2,3)*HH(3,2) )
        rH(1,2)=rdet*( HH(3,2)*HH(1,3)-HH(3,3)*HH(1,2) )
        rH(1,3)=rdet*( HH(1,2)*HH(2,3)-HH(1,3)*HH(2,2) )
        rH(2,1)=rdet*( HH(2,3)*HH(3,1)-HH(2,1)*HH(3,3) )
        rH(2,2)=rdet*( HH(3,3)*HH(1,1)-HH(3,1)*HH(1,3) )
        rH(2,3)=rdet*( HH(1,3)*HH(2,1)-HH(1,1)*HH(2,3) )
        rH(3,1)=rdet*( HH(2,1)*HH(3,2)-HH(2,2)*HH(3,1) )
        rH(3,2)=rdet*( HH(3,1)*HH(1,2)-HH(3,2)*HH(1,1) )
        rH(3,3)=rdet*( HH(1,1)*HH(2,2)-HH(1,2)*HH(2,1) )

        return
        end

!########################################################################
        subroutine creatematrix(  &
                           box,alpha,beta,gamma,HH,rH,tHHH)
!########################################################################
        implicit none
        real(8) :: box(3)
        real(8) :: alpha,beta,gamma
        real(8) :: HH(3,3),rH(3,3)
	real(8) :: tH(3,3), tHHH(3,3)
 	real(8) :: av(3),bv(3),cv(3)
	real(8) :: wa
	integer(4) :: i,j,k,l,m,n

!### create HH, rHH, tHHH matrix ###

        av(1)=box(1)
        av(2)=0d0
        av(3)=0d0
        bv(1)=box(2)*dcos(gamma)
        bv(2)=box(2)*dsin(gamma)
        bv(3)=0d0
        cv(1)=box(3)*dcos(beta)
        cv(2)=box(3)/dsin(gamma)*(                                &
               dcos(alpha)*dsin(gamma)**2                         &
              -dcos(gamma)*(dcos(beta)-dcos(alpha)*dcos(gamma))   &
              )
        cv(3)=box(3)*dsqrt(                                       & 
               dsin(alpha)**2 - (dcos(beta)-dcos(alpha)*dcos(gamma) )**2  &
               /dsin(gamma)**2 )
        HH(1,1)=av(1)
        HH(2,1)=av(2) ! 0d0
        HH(3,1)=av(3) ! 0d0
        HH(1,2)=bv(1)
        HH(2,2)=bv(2)
        HH(3,2)=bv(3) ! 0d0
        HH(1,3)=cv(1)
        HH(2,3)=cv(2)
        HH(3,3)=cv(3)
        call reverse(HH,rH)

!### transverse of HH ###
	do i=1,3
	do j=1,3
	  tH(i,j)=HH(j,i)
	enddo
	enddo

        tHHH=0d0
	do i=1,3
	do j=1,3
	  wa=0d0
	  do k=1,3
	    wa=wa+tH(i,k)*HH(k,j)
	  enddo
	  tHHH(i,j)=wa
	enddo
	enddo

        return
        end

