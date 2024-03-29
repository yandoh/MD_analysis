!Copyright (c) 2022 Yoshimichi ANDOH
!Released under the MIT license
!https://opensource.org/licenses/mit-license.php

! === inputs ===
! (1) ./sys_info
! (2) ./massinfo.mdff
! (3) ./DCD
!
      implicit none
!
      integer(4),allocatable::nav(:),nmv(:)
      integer(4)::h,i,k,m,komtot,kom
      integer(4)::komi,komj,hi,hj
      integer(4)::iser
      real(8)::dx,dy,dz,r2,Rg,Ri
      real(8),allocatable::sRg(:),sRi(:)
      real(8),allocatable::sRg2(:),sRi2(:)
      real(8),allocatable::sdvg(:),sdvi(:)
      integer(4) :: totnp
      real(8) :: alpha,beta,gamma,box(3)
      integer(4) :: molcount
      real(8),allocatable :: mass(:,:), molmass(:)
      real(8),allocatable :: pos(:,:,:)
      real(8),allocatable :: cntmol0(:,:,:)
      real(8) :: pi
!dcd
      character(len=4)::Aname
      integer(4) :: nstr, nptot, iflame, nflame
      integer(4),dimension(20)::icntrl
      real(8),dimension(6)::cellstr
      real(4),allocatable::flpx(:),flpy(:),flpz(:)
!dcd

	molcount=0
	pi=dacos(-1d0)
	totnp=0
!### input sysinfo ###
	open(21,file='sys_info',status='old')
        read(21,*) ! skip
 	read(21,*) komtot
 	allocate( nav(komtot), nmv(komtot) )
        allocate( sRg(komtot), sRi(komtot) )
        allocate( sRg2(komtot), sRi2(komtot) )
        allocate( sdvg(komtot), sdvi(komtot) )
        do kom=1,komtot
          read(21,*) nav(kom),nmv(kom)
          molcount=molcount+nmv(kom)
	  totnp =totnp +nmv(kom)*nav(kom)
        enddo
	close(21)

        sRg=0d0; sRg2=0d0
        sRi=0d0; sRi2=0d0

!allocate arrays
	allocate(mass(komtot,totnp), molmass(komtot))
	allocate(pos(3,komtot,totnp))
        allocate(cntmol0(3,komtot,molcount))
	
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
! by catdcd
    if(abs(alpha).le.1d0)then  ! cos(angle)
      alpha=acos(alpha)   ! rad
      beta =acos(beta)    ! rad
      gamma=acos(gamma)   ! rad
    else  !! [degree]
      alpha=alpha/180d0*pi
      beta=beta/180d0*pi
      gamma=gamma/180d0*pi
    endif

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

!	^^^ calc Rg, Ri ^^^!
	DO KOM=1,komtot
	do h=1,nmv(KOM)
	  i=(h-1)*nav(KOM)
          Rg=0d0
          Ri=0d0
	  do m=1,nav(KOM)
            dx=pos(1,kom,i+m)-cntmol0(1,kom,h)
            dy=pos(2,kom,i+m)-cntmol0(2,kom,h)
            dz=pos(3,kom,i+m)-cntmol0(3,kom,h)
            r2=dx**2+dy**2+dz**2
            Rg=Rg+mass(kom,m)*r2
            Ri=Ri+r2
	  enddo
          Rg=Rg/molmass(kom)
          Ri=Ri/nav(KOM)
          sRg(kom)=sRg(kom)+sqrt(Rg)
          sRi(kom)=sRi(kom)+sqrt(Ri)
!!        write(*,*) sqrt(Rg) !; stop
          sRg2(kom)=sRg2(kom)+Rg
          sRi2(kom)=sRi2(kom)+Ri
	enddo
	ENDDO
        write(10,*) iflame, '-th flame analyzed/', nflame
        call flush(10)
      ENDDO ! iflame
!stop

      close(34)
      close(10)

!### average over flame ###!
      do KOM=1,komtot
!!      write(*,*) sRg(kom),nmv(kom)*nflame
        sRg(kom)=sRg(kom)/nmv(kom)/nflame
        sRi(kom)=sRi(kom)/nmv(kom)/nflame
        sRg2(kom)=sRg2(kom)/nmv(kom)/nflame
        sRi2(kom)=sRi2(kom)/nmv(kom)/nflame
        sdvg(kom)=sqrt(sRg2(kom)-sRg(kom)**2)
        sdvi(kom)=sqrt(sRi2(kom)-sRi(kom)**2)
      enddo

!###  output  ###
      write(*,'(a1,i5,i10)') '#', komtot, nflame
      do kom=1,komtot
	write(*,'(i5,2f10.3,a,2f10.3,a)') &
    &               kom, sRg(kom), sRi(kom), &
    &              " (", sdvg(kom), sdvi(kom), ")"
      enddo

      stop
      end
