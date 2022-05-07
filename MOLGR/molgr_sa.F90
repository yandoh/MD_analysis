!Copyright (c) 2022 Yoshimichi ANDOH
!Released under the MIT license
!https://opensource.org/licenses/mit-license.php

! Stand alone version of molgr.F90
!
! === inputs ===
! (1) ./sys_info
! (2) ./massinfo.mdff
! (3) ./DCD
!
      implicit none
      integer(4),parameter::ncopy=1
!
      integer(4),allocatable::nav(:),nmv(:)
      integer(4)::h,i,k,m,komtot,kom
      integer(4)::komi,komj,hi,hj
      integer(4)::iser
      integer(4),parameter :: isermax=(ncopy+1+ncopy)**3
      real(8)::xij,yij,zij,rij,rij2
      integer(4) :: totnp
      real(8) :: alpha,beta,gamma,box(3)
      real(8),allocatable :: cntmol(:,:,:,:)
      integer(4) :: isermin0
      real(8),allocatable :: hst(:,:),gr(:,:)
      real(8),allocatable :: rho(:)
      real(8) :: dr,rmax,rmax2
      integer(4) :: ir,ndr
      integer(4) :: molcount
      real(8) :: radius,dvolume,volume,svolume
      real(8) :: pi
      real(8),allocatable :: mass(:,:), molmass(:)
      integer(4) ix,iy,iz
      real(8)::avec(3),bvec(3),cvec(3)
!     real(8) :: HH(3,3), rH(3,3), tHHH(3,3)
      real(8),allocatable :: pos(:,:,:)
      real(8),allocatable :: cntmol0(:,:,:)
!dcd
      character(len=4)::Aname
      integer(4) :: nstr, nptot, iflame, nflame
      integer(4),dimension(20)::icntrl
      real(8),dimension(6)::cellstr
      real(4),allocatable::flpx(:),flpy(:),flpz(:)
!dcd

!parameter for g(r)
	rmax=25d0          !! max of r value
	rmax2=rmax*rmax
	dr=0.1d0           !! delta r
	ndr=int(rmax/dr)
	pi=dacos(-1d0)

	molcount=0
	volume=0d0
	totnp=0
!### input sysinfo ###
	open(21,file='sys_info',status='old')
        read(21,*) ! skip
 	read(21,*) komtot
 	allocate( nav(komtot), nmv(komtot) )
        do kom=1,komtot
          read(21,*) nav(kom),nmv(kom)
          molcount=molcount+nmv(kom)
	  totnp =totnp +nmv(kom)*nav(kom)
        enddo
	close(21)

	isermin0=(ncopy+1+ncopy)**3/2+1

!allocate arrays
	allocate( hst(ndr,komtot) )
	allocate( gr(ndr,komtot) )
	allocate( rho(komtot) )
	allocate( mass(komtot,totnp), molmass(komtot) )
	allocate( pos(3,komtot,totnp) )
        allocate( cntmol0(3,komtot,molcount) )
        allocate( cntmol(3,komtot,molcount,(ncopy+1+ncopy)**3) )
	hst=0d0
	gr=0d0
	
!### input mass ###
	open(3,file='massinfo.mdff',status='old')
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

!### output log ###
      open(10,file='log',status='replace')
      write(10,*) komtot
      write(10,*) molmass

!### input .dcd file & analyze hst(r) ###
      open(34,file='DCD',form='unformatted',status='old')
      read(34) Aname,icntrl
      read(34) nstr
      read(34) nptot
      allocate( flpx(nptot),flpy(nptot),flpz(nptot) )

      nflame=icntrl(1)
      DO iflame=1,nflame

	read(34) cellstr
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
      gamma =cellstr(2)   ! cos(gamma)   by catdcd
      box(2)=cellstr(3)   ! |b|
      beta  =cellstr(4)   ! cos(beta)    by catdcd
      alpha =cellstr(5)   ! cos(alpha)   by catdcd
      box(3)=cellstr(6)   ! |c|
! by DCDMOL
      alpha=alpha/180d0*pi
      beta =beta /180d0*pi
      gamma=gamma/180d0*pi
! by catdcd
!     alpha=acos(alpha)   ! rad
!     beta =acos(beta)    ! rad
!     gamma=acos(gamma)   ! rad

      call createvectors(box,alpha,beta,gamma,avec,bvec,cvec)

!### calc. COM ###
	cntmol0=0d0
	DO KOM=1,komtot
	do h=1,nmv(KOM)
	i=(h-1)*nav(KOM)
	do k=1,3  ! x, y, z
	  do m=1,nav(KOM)
	    cntmol0(k,kom,h)=cntmol0(k,kom,h)+pos(k,kom,i+m)*mass(kom,m)
	  enddo
	  cntmol0(k,kom,h)=cntmol0(k,kom,h)/molmass(kom)
	enddo
	enddo
	ENDDO

!	^^^ record cntmol ^^^
        iser=0
 	do kom=1,komtot
 	do ix=-ncopy,+ncopy ; do iy=-ncopy,+ncopy ; do iz=-ncopy,+ncopy
        iser=iser+1
 	do h=1,nmv(kom)
	  cntmol(:,kom,h,iser)=cntmol0(:,kom,h)  &
     &                        +ix*avec(:) + iy*bvec(:) + iz*cvec(:)
 	enddo ! h
	enddo; enddo; enddo
 	enddo ! kom

        call cellvolume(avec,bvec,cvec,volume)
	svolume=svolume+volume

!	^^^ calc hst(r) ^^^!
 	do komi=1,komtot
 	do hi=1,nmv(komi)
	do iser=1,isermax
 	do komj=1,komtot
 	do hj=1,nmv(komj)
 	    xij=cntmol(1,komj,hj,iser)-cntmol(1,komi,hi,isermin0)
 	    yij=cntmol(2,komj,hj,iser)-cntmol(2,komi,hi,isermin0)
   	    zij=cntmol(3,komj,hj,iser)-cntmol(3,komi,hi,isermin0)
	    rij2=xij**2+yij**2+zij**2
	  if(rij2==0d0) cycle
	  if(rij2.gt.rmax2) cycle
	    rij=sqrt(rij2)
	    ir=dint(rij/dr)+1
	    hst(ir,komi)=hst(ir,komi)+1d0
	enddo ! h
	enddo ! kom
	enddo ! iser
	enddo ! h
	enddo ! kom

        write(10,*) iflame, '-th flame analyzed/', nflame
        call flush(10)
      ENDDO ! iflame

      close(34)
      close(10)

!### average over flame ###!
      hst=hst/dble(nflame)
      svolume=svolume/dble(nflame)

!### calc gr ###!
      do kom=1,komtot
        rho(kom)=nmv(kom)/svolume 
      enddo

      do kom=1,komtot
      do ir=1,ndr
 	radius=dr*(ir-0.5d0)
	dvolume=(4d0*pi*radius**2)*dr
 	gr(ir,kom)=hst(ir,kom)/nmv(kom)/(dvolume*rho(kom))  
      enddo
      enddo

!###  output  ###
      write(*,'(a1,e22.15,2i5,i10)') '#', svolume, komtot, ndr, nflame
      write(*,'(a1,99e22.15)') '#', rho(1:komtot)
      do ir=1,ndr
 	radius=dr*(ir-0.5d0)
	write(*,'(f12.5,99e20.10)') radius,gr(ir,1:komtot)
      enddo

      stop
      end

!############################################################
        subroutine cellvolume(av,bv,cv,volume)
!########################################################################
        implicit none
        real(8)::av(3),bv(3),cv(3),volume
        real(8)::bvcv(3), avdotbvcv
        integer(4)::k

        bvcv(1)=bv(2)*cv(3)-bv(3)*cv(2)
        bvcv(2)=bv(3)*cv(1)-bv(1)*cv(3)
        bvcv(3)=bv(1)*cv(2)-bv(2)*cv(1)

        avdotbvcv=0d0
        do k=1,3
          avdotbvcv=avdotbvcv+av(k)*bvcv(k)
        enddo

        volume=avdotbvcv

        endsubroutine cellvolume

!########################################################################
        subroutine createvectors(  &
                           box,alpha,beta,gamma,av,bv,cv)
!########################################################################
        implicit none
        real(8) :: box(3)
        real(8) :: alpha,beta,gamma
        real(8) :: av(3),bv(3),cv(3)

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

        return
        end subroutine createvectors
