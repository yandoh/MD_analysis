!Copyright (c) 2022 Yoshimichi ANDOH
!Released under the MIT license
!https://opensource.org/licenses/mit-license.php

! Stand alone version of partgr.F90
!
! === inputs ===
! (1) ./sys_info
! (2) ./massinfo.mdff
! (3) ./molparts
! (4) ./DCD
!
      implicit none
      integer(4),parameter :: ncopy=1
      integer(4),parameter :: isermax=(ncopy+1+ncopy)**3
      integer(4),parameter :: maxjatoms=100 ! maxium number of atoms in a part of molecule
!
      integer(4),allocatable :: nav(:),nmv(:)
      integer(4),allocatable :: natoms_in_part(:)
      integer(4),allocatable :: kom_part(:)
      integer(4),allocatable :: atom_members(:,:)
      real(8),allocatable :: cntpart(:,:,:,:)
      real(8),allocatable :: hst(:,:),gr(:,:)
      real(8),allocatable :: rho(:)
      real(8),allocatable :: mass(:,:)
      real(8),allocatable :: pos(:,:,:)
      real(8),allocatable :: cntpart0(:,:,:)
      real(8),allocatable :: partmass(:)
!
      integer(4) h,i,j,k,L,m,komtot,kom
      integer(4) hi,hj,L1,L2,kom1,kom2
      integer(4) iser
      integer(4) npmax
      integer(4) isermin0
      integer(4) ir,ndr
      integer(4) molmax
      integer(4) ix,iy,iz
      integer(4) nparts,id,idparts
      real(8) xij,yij,zij,rij,rij2
      real(8) alpha,beta,gamma,box(3)
      real(8) dr,rmax,rmax2
      real(8) radius,dvolume,volume,svolume
      real(8) pi
      real(8) avec(3),bvec(3),cvec(3)
!dcd
      character(len=4)::Aname
      integer(4) :: nstr, nptot, iflame, nflame
      integer(4),dimension(20)::icntrl
      real(8),dimension(6)::cellstr
      real(4),allocatable::flpx(:),flpy(:),flpz(:)
!dcd

!parameter for g(r)
	rmax=40d0          !! max of r value
	rmax2=rmax*rmax
	dr=0.1d0           !! delta r
	ndr=int(rmax/dr)
	pi=dacos(-1d0)

	molmax=0
	svolume=0d0
	npmax=0
!### input sysinfo ###
	open(21,file='sys_info',status='old')
        read(21,*) ! skip
 	read(21,*) komtot
 	allocate( nav(komtot), nmv(komtot) )
        do kom=1,komtot
          read(21,*) nav(kom),nmv(kom)
          molmax=max(molmax,nmv(kom))
	  npmax =max(npmax,nmv(kom)*nav(kom))
        enddo
	close(21)

	isermin0=(ncopy+1+ncopy)**3/2+1

!allocate arrays
	allocate( mass(komtot,npmax) )
	allocate( pos(3,komtot,npmax) )
	
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

!### input parts ###
	open(14,file='molparts',status='old')
        read(14,*) nparts
        allocate( natoms_in_part(nparts) )
        allocate( kom_part(nparts) )
        allocate( atom_members(maxjatoms,nparts) )
        do i=1,nparts 
          read(14,*) kom_part(i), natoms_in_part(i), &
          (atom_members(j,i),j=1,natoms_in_part(i))
        enddo
        close(14)

        allocate( partmass(nparts) )
        allocate( cntpart0(3,molmax,nparts) )
        allocate( cntpart(3,molmax,nparts,(ncopy+1+ncopy)**3) )

        idparts=0
        do i=1,nparts
        do j=i,nparts
          idparts=idparts+1
        enddo
        enddo
        allocate( rho(idparts) )
	allocate( hst(ndr,idparts) )
	allocate( gr(ndr,idparts) )
	hst=0d0
	gr=0d0

!### calc parts mass ###
        DO L=1,nparts
          KOM=kom_part(L)+1
          partmass(L)=0d0
          do j=1,natoms_in_part(L)
            m=atom_members(j,L)+1
            partmass(L)=partmass(L)+mass(KOM,m)
	  enddo
        ENDDO

!### output log ###
      open(10,file='log',status='replace')
      write(10,*) komtot
      write(10,*) nparts,idparts
      write(10,*) kom_part
      write(10,*) natoms_in_part

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

!### calc. COM of each part ###
	cntpart0=0d0
        DO L=1,nparts
          KOM=kom_part(L)+1
	  do h=1,nmv(KOM)
	    i=(h-1)*nav(KOM)
            do j=1,natoms_in_part(L)
              m=atom_members(j,L)+1
              cntpart0(:,h,L)=cntpart0(:,h,L) &
     &                         +pos(:,kom,i+m)*mass(kom,m)
	    enddo
	    cntpart0(:,h,L)=cntpart0(:,h,L)/partmass(L)
          enddo
	ENDDO

!	^^^ record cntpart ^^^
        do L=1,nparts
        KOM=kom_part(L)+1
        iser=0
 	do ix=-ncopy,+ncopy ; do iy=-ncopy,+ncopy ; do iz=-ncopy,+ncopy
        iser=iser+1
 	do h=1,nmv(kom)
	  cntpart(:,h,L,iser)=cntpart0(:,h,L)  &
     &                        +ix*avec(:) + iy*bvec(:) + iz*cvec(:)
 	enddo ! h
	enddo; enddo; enddo
 	enddo ! kom

        call cellvolume(avec,bvec,cvec,volume)
	svolume=svolume+volume

!	^^^ calc hst(r) ^^^!
        id=0
        DO L1=1,nparts
        DO L2=L1,nparts
        id=id+1
        kom1=kom_part(L1)+1
        kom2=kom_part(L2)+1
 	do hi=1,nmv(kom1)
	do iser=1,isermax
 	do hj=1,nmv(kom2)
 	    xij=cntpart(1,hj,L2,iser)-cntpart(1,hi,L1,isermin0)
 	    yij=cntpart(2,hj,L2,iser)-cntpart(2,hi,L1,isermin0)
   	    zij=cntpart(3,hj,L2,iser)-cntpart(3,hi,L1,isermin0)
	    rij2=xij**2+yij**2+zij**2
	  if(rij2==0d0) cycle
	  if(rij2.gt.rmax2) cycle
	    rij=sqrt(rij2)
	    ir=dint(rij/dr)+1
	    hst(ir,id)=hst(ir,id)+1d0
	enddo ! hj
	enddo ! iser
	enddo ! hi
        enddo ! L2
        enddo ! L1

        write(10,*) iflame, '-th flame analyzed/', nflame
        call flush(10)
      ENDDO ! iflame

      close(34)
      close(10)

!### average over flame ###!
      hst=hst/dble(nflame)
      svolume=svolume/dble(nflame)

!### calc gr ###!
      id=0
      do L1=1,nparts
      do L2=L1,nparts
        id=id+1
        kom1=kom_part(L1)+1
        kom2=kom_part(L2)+1
        rho(id)=nmv(kom2)/svolume 
        do ir=1,ndr
 	  radius=dr*(ir-0.5d0)
	  dvolume=(4d0*pi*radius**2)*dr
 	  gr(ir,id)=hst(ir,id)/nmv(kom1)/(dvolume*rho(id))
        enddo !ir
      enddo ! L2
      enddo ! L1

!###  output  ###
      write(*,'(a1,e22.15,2i5,i10)') '#', svolume, idparts, ndr, nflame
      write(*,'(a1,99e22.15)') '#', rho
      do ir=1,ndr
 	radius=dr*(ir-0.5d0)
	write(*,'(f12.5,99e20.10)') radius,(gr(ir,L),L=1,idparts)
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
