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
      real(8),parameter :: dt=1d-15
      integer(4),parameter :: dcd_interval=500
      integer(4),parameter::ncopy=1
!
      integer(4),allocatable::nav(:),nmv(:)
      integer(4)::h,i,k,m,komtot,kom
      integer(4)::iser
      integer(4),parameter :: isermax=(ncopy+1+ncopy)**3
      integer(4) :: totnp
      real(8) :: alpha,beta,gamma,box(3)
      real(8),allocatable :: cntmol(:,:,:,:,:)
      integer(4) :: isermin0
      real(8) :: dr2,dr2tmp
      integer(4) :: molcount
      real(8),allocatable :: mass(:,:), molmass(:)
      integer(4) ix,iy,iz
      real(8)::avec(3),bvec(3),cvec(3)
      real(8)::dx,dy,dz,time
      real(8),allocatable :: pos(:,:,:)
      real(8),allocatable :: cntmol0(:,:,:)
      real(8) :: pi
!dcd
      character(len=4)::Aname
      integer(4) :: nstr, nptot, iflame, jflame, dflame, nflame
      integer(4),dimension(20)::icntrl
      real(8),dimension(6)::cellstr
      real(4),allocatable::flpx(:),flpy(:),flpz(:)
!dcd
      real(8),allocatable::msd(:,:)
      real(8),allocatable::dmsd(:)
      integer(4),allocatable::icount(:,:)

	pi=dacos(-1d0)
	molcount=0
	totnp=0
!### input sysinfo ###
	open(21,file='sys_info',status='old')
        read(21,*) ! skip
 	read(21,*) komtot
 	allocate( nav(komtot), nmv(komtot) )
        do kom=1,komtot
          read(21,*) nav(kom),nmv(kom)
          molcount=max(molcount,nmv(kom))
	  totnp =max(totnp,nmv(kom)*nav(kom))
        enddo
	close(21)

	isermin0=(ncopy+1+ncopy)**3/2+1

!allocate arrays
	allocate( mass(komtot,totnp), molmass(komtot) )
	allocate( pos(3,komtot,totnp) )
        allocate( cntmol0(3,molcount,komtot) )
	
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

!### input .dcd file ###
      open(34,file='DCD',form='unformatted',status='old')
      read(34) Aname,icntrl
      read(34) nstr
      read(34) nptot
      allocate( flpx(nptot),flpy(nptot),flpz(nptot) )
      nflame=icntrl(1)
      allocate( cntmol(3,molcount,komtot,isermax,nflame) )
      allocate( msd(0:nflame,komtot), icount(0:nflame,komtot) )
      allocate( dmsd(komtot) )
      msd=0d0
      icount=0

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

!### calc. COM ###
	cntmol0=0d0
	DO KOM=1,komtot
	do h=1,nmv(KOM)
	i=(h-1)*nav(KOM)
	do k=1,3  ! x, y, z
	  do m=1,nav(KOM)
	    cntmol0(k,h,kom)=cntmol0(k,h,kom)+pos(k,kom,i+m)*mass(kom,m)
	  enddo
          cntmol0(k,h,kom)=cntmol0(k,h,kom)/molmass(kom)
	enddo
	enddo
	ENDDO

!	^^^ record cntmol ^^^
        iser=0
 	do kom=1,komtot
 	do ix=-ncopy,+ncopy ; do iy=-ncopy,+ncopy ; do iz=-ncopy,+ncopy
        iser=iser+1
 	do h=1,nmv(kom)
          cntmol(:,h,kom,iser,iflame)=cntmol0(:,h,kom)  &
     &                        +ix*avec(:) + iy*bvec(:) + iz*cvec(:)
 	enddo ! h
        enddo; enddo; enddo
 	enddo ! kom

      ENDDO ! iflame

      close(34)

!### take auto-correlation ###!
      do iflame=1,nflame-1
      do jflame=iflame,nflame
	DO KOM=1,komtot
	do h=1,nmv(KOM)
          dr2=1d+10  ! very large inital value
        do iser=1,isermax
          dx=cntmol(1,h,kom,iser,jflame)-cntmol(1,h,kom,isermin0,iflame)
          dy=cntmol(2,h,kom,iser,jflame)-cntmol(2,h,kom,isermin0,iflame)
          dz=cntmol(3,h,kom,iser,jflame)-cntmol(3,h,kom,isermin0,iflame)
          dr2tmp=dx**2+dy**2+dz**2
          if(dr2tmp .lt. dr2)then
            dr2=dr2tmp  !! select minimum dr
          endif
        enddo
        dflame=jflame-iflame
        msd(dflame,kom)=msd(dflame,kom)+dr2
        icount(dflame,kom)=icount(dflame,kom)+1
        enddo !h 
        ENDDO !KOM
      enddo
      enddo

!### average msd ###!
      do KOM=1,komtot
      do dflame=0,nflame-1
        if(icount(dflame,kom) /= 0)then
          msd(dflame,kom)=msd(dflame,kom)/dble(icount(dflame,kom))
        endif
      enddo
      enddo

!### output msd ###
      write(*,'(a,es10.3)') '# parameter: dt          =', dt
      write(*,'(a,i10)')    '# parameter: dcd_interval=', dcd_interval
      write(*,'(a,i5,i10)') '#', komtot, nflame
      write(*,'(a       )') '#   t [ps]    D [cm^2/s]    msd   [A^2]'
      do dflame=0,nflame-1
        time=dt*dcd_interval*dflame*1d+12
        if(dflame==0)then
          dmsd=0d0
        else
          dmsd(:)=msd(dflame,:)-msd(dflame-1,:)
        endif
        dmsd=dmsd*1d-20/(dt*dcd_interval)/6d0 ! [m^2/s]
        dmsd=dmsd*1d+4 ! [cm^2/s]   1m=100cm
        write(*,'(f9.3,99es15.7)') time, dmsd(:), msd(dflame,:)
      enddo

      close(10)

      stop
      end

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
