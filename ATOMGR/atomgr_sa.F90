!Copyright (c) 2022 Yoshimichi ANDOH
!Released under the MIT license
!https://opensource.org/licenses/mit-license.php

! Stand alone version of atomgr.F90
!
! === inputs ===
! (1) ./sys_info
! (2) ./atompairs
! (3) ./DCD
!
      implicit none
      integer(4),parameter::ncopy=1
!
      integer(4) :: natompairs
      integer(4),allocatable :: pairid(:,:)
!
      integer(4),allocatable::nav(:),nmv(:)
      integer(4)::h,i,j,k,L,m,komtot,kom
      integer(4)::komi,komj,hi,hj,atmi,atmj
      integer(4)::iser
      integer(4),parameter :: isermax=(ncopy+1+ncopy)**3
      real(8)::xij,yij,zij,rij,rij2
      integer(4) :: totnp
      real(8) :: alpha,beta,gamma,box(3)
      integer(4) :: isermin0
      real(8),allocatable :: hst(:,:),gr(:,:)
      real(8) :: dr,rmax,rmax2
      integer(4) :: ir,ndr
      real(8),allocatable :: rho(:)
      real(8) :: radius,dvolume,volume,svolume
      real(8) :: pi
      integer(4) ix,iy,iz
      real(8)::avec(3),bvec(3),cvec(3)
      real(8),allocatable :: pos0(:,:,:)
      real(8),allocatable :: pos(:,:,:,:)
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

	volume=0d0
	totnp=0
!### input sysinfo ###
	open(21,file='sys_info',status='old')
        read(21,*) ! skip
 	read(21,*) komtot
 	allocate( nav(komtot), nmv(komtot) )
        do kom=1,komtot
          read(21,*) nav(kom),nmv(kom)
	  totnp =totnp +nmv(kom)*nav(kom)
        enddo
	close(21)

	isermin0=(ncopy+1+ncopy)**3/2+1

!### input atompairs ###
	open(22,file='atompairs',status='old')
        read(22,*) natompairs
        allocate( pairid(4,natompairs) )
 	allocate( rho(natompairs) )
        do i=1,natompairs
          read(22,*) pairid(1:4,i)
	enddo

!allocate arrays
	allocate( hst(ndr,natompairs) )
	allocate( gr(ndr,natompairs) )
	allocate( pos0(3,totnp,komtot) )
	allocate( pos(3,totnp,komtot,(ncopy+1+ncopy)**3) )
	hst=0d0
	gr=0d0
	
!### output log ###
      open(10,file='log',status='replace')
      write(10,*) komtot
      write(10,*) natompairs

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
            pos0(1,i+m,kom)=flpx(k)
            pos0(2,i+m,kom)=flpy(k)
            pos0(3,i+m,kom)=flpz(k)
          enddo
        enddo
	enddo

      box(1)=cellstr(1)   ! |a|
      gamma =cellstr(2)   ! cos(gamma)   by catdcd
      box(2)=cellstr(3)   ! |b|
      beta  =cellstr(4)   ! cos(beta)    by catdcd
      alpha =cellstr(5)   ! cos(alpha)   by catdcd
      box(3)=cellstr(6)   ! |c|

    if(abs(alpha).le.1d0)then  ! cos(angle)
      alpha=acos(alpha)   ! rad
      beta =acos(beta)    ! rad
      gamma=acos(gamma)   ! rad
    else  !! [degree]
      alpha=alpha/180d0*pi
      beta=beta/180d0*pi
      gamma=gamma/180d0*pi
    endif

      call createvectors(box,alpha,beta,gamma,avec,bvec,cvec)
      call cellvolume(avec,bvec,cvec,volume)
      svolume=svolume+volume

!	^^^ record pos ^^^
        iser=0
 	do kom=1,komtot
 	do ix=-ncopy,+ncopy ; do iy=-ncopy,+ncopy ; do iz=-ncopy,+ncopy
        iser=iser+1
 	do h=1,nmv(kom)
          i=(h-1)*nav(kom)
          do m=1,nav(kom)
	    pos(:,i+m,kom,iser)=pos0(:,i+m,kom)  &
     &                         +ix*avec(:) + iy*bvec(:) + iz*cvec(:)
 	enddo ! m
 	enddo ! h
	enddo; enddo; enddo
 	enddo ! kom

      DO L=1,natompairs
        komi=pairid(1,L)+1  ! id start from 0  in atompairs file
        komj=pairid(2,L)+1  ! id start from 0  in atompairs file
        atmi=pairid(3,L)+1  ! id start from 0  in atompairs file
        atmj=pairid(4,L)+1  ! id start from 0  in atompairs file
!	^^^ calc hst(r) ^^^!
 	do hi=1,nmv(komi)
          i=(hi-1)*nav(komi)+atmi
	do iser=1,isermax
 	do hj=1,nmv(komj)
          j=(hj-1)*nav(komj)+atmj
          if(komi==komj .and. hi==hj) cycle
 	    xij=pos(1,j,komj,iser)-pos(1,i,komi,isermin0)
 	    yij=pos(2,j,komj,iser)-pos(2,i,komi,isermin0)
   	    zij=pos(3,j,komj,iser)-pos(3,i,komi,isermin0)
	    rij2=xij**2+yij**2+zij**2
	  if(rij2==0d0) cycle
	  if(rij2.gt.rmax2) cycle
	    rij=sqrt(rij2)
	    ir=dint(rij/dr)+1
	    hst(ir,L)=hst(ir,L)+1d0
	enddo ! hj
	enddo ! iser
	enddo ! hi
      ENDDO

        write(10,*) iflame, '-th flame analyzed/', nflame
        call flush(10)
      ENDDO ! iflame

      close(34)
      close(10)

!### average over flame ###!
      hst=hst/dble(nflame)
      svolume=svolume/dble(nflame)

!### calc gr ###!
      do L=1,natompairs
        komj=pairid(2,L)+1
        rho(L)=nmv(komj)/svolume 
      enddo

      do L=1,natompairs
      komi=pairid(1,L)+1
      komj=pairid(2,L)+1
      do ir=1,ndr
 	radius=dr*(ir-0.5d0)
	dvolume=(4d0*pi*radius**2)*dr
 	gr(ir,L)=hst(ir,L)/nmv(komi)/(dvolume*rho(L))  
      enddo ! ir
      enddo ! L

!###  output  ###
      write(*,'(a1,e22.15,2i5,i10)') '#', &
     &                  svolume, natompairs, ndr, nflame
      write(*,'(a1,99e22.15)') '#', rho
      do ir=1,ndr
 	radius=dr*(ir-0.5d0)
 	write(*,'(f12.5,99e20.10)') radius,gr(ir,:)
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
