!Copyright (c) 2022 Yoshimichi ANDOH
!Released under the MIT license
!https://opensource.org/licenses/mit-license.php

!
! === inputs ===
! (1) ./sys_info
! (2) ./massinfo.mdff
! (3) ./vectors
! (4) ./molparts
! (5) ./DCD (wrapped by DCDMOL)
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
      real(8),allocatable :: mass(:,:)
      real(8),allocatable :: pos(:,:,:)
      real(8),allocatable :: dpos(:,:,:)
      real(8),allocatable :: cntpart0(:,:,:)
      real(8),allocatable :: partmass(:)
!
      integer(4) h,i,j,k,L,m,komtot,kom
      integer(4) hi,hj
      integer(4) iser
      integer(4) npmax
      integer(4) isermin0
      integer(4) molmax
      integer(4) ix,iy,iz
      integer(4) nparts
      real(8) xij,yij,zij,rij,rij2
      real(8) alpha,beta,gamma,box(3)
      real(8) pi
      real(8) avec(3),bvec(3),cvec(3)
      integer(4) nvectors
      integer(4) ia(2),ja(2)
      real(8) Rmin, Rmin2
      integer(4) ii1,ij1,ii2,ij2,ji1,jj1,ji2,jj2
      real(8) veci1(3),veci2(3),veci3(3)
      real(8) vecj1(3),vecj2(3),vecj3(3)
      real(8) costheta,theta, absvi, absvj, thetarad, dthetarad, prbs
      integer(4),allocatable :: hsttheta(:)
      real(8),allocatable :: prbtheta(:)
      real(8),parameter :: dtheta=5d0 ! degree
      integer(4) ntheta,itheta,nsample
      real(8) :: Sv=0d0
!dcd
      character(len=4)::Aname
      integer(4) :: nstr, nptot, iflame, nflame
      integer(4),dimension(20)::icntrl
      real(8),dimension(6)::cellstr
      real(4),allocatable::flpx(:),flpy(:),flpz(:)
!dcd

!### set constants ###
	pi=dacos(-1d0)
        ntheta=int(180d0/dtheta)
        allocate(hsttheta(1:ntheta))
        allocate(prbtheta(1:ntheta))
        hsttheta=0
        prbtheta=0d0

	molmax=0
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
	allocate( pos(3,npmax,komtot) )
	allocate( dpos(3,npmax,komtot) )
	
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

!### input vectors ###
	open(15,file='vectors',status='old')
        read(15,*) nvectors, Rmin
        read(15,*) ia(1), ja(1)
        read(15,*) ia(2), ja(2)
        close(15)

        Rmin2=Rmin**2

!### calc parts mass ###
        DO L=1,nparts
          KOM=kom_part(L)+1
          partmass(L)=0d0
          do j=1,natoms_in_part(L)
            m=atom_members(j,L) +1
            partmass(L)=partmass(L)+mass(KOM,m)
	  enddo
        ENDDO

!### output log ###
      open(10,file='log',status='replace')
      write(10,*) komtot
      write(10,*) nparts
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
            pos(1,i+m,kom)=flpx(k)
            pos(2,i+m,kom)=flpy(k)
            pos(3,i+m,kom)=flpz(k)
          enddo
        enddo
	enddo

      box(1)=cellstr(1)   ! |a|
      gamma =cellstr(2)   ! cos(gamma)   by catdcd
      box(2)=cellstr(3)   ! |b|
      beta  =cellstr(4)   ! cos(beta)    by catdcd
      alpha =cellstr(5)   ! cos(alpha)   by catdcd
      box(3)=cellstr(6)   ! |c|

      alpha=alpha/180d0*pi
      beta =beta /180d0*pi
      gamma=gamma/180d0*pi
!      alpha=acos(alpha)   ! rad
!      beta =acos(beta)    ! rad
!      gamma=acos(gamma)   ! rad

      call createvectors(box,alpha,beta,gamma,avec,bvec,cvec)

!### calc. COM of each part ###
	cntpart0=0d0
        DO L=1,nparts
          KOM=kom_part(L)+1
	  do h=1,nmv(KOM)
	    i=(h-1)*nav(KOM)
            do j=1,natoms_in_part(L)
              m=atom_members(j,L) +1
              cntpart0(:,h,L)=cntpart0(:,h,L) &
     &                         +pos(:,i+m,kom)*mass(kom,m)
	    enddo
!calc CoM
	    cntpart0(:,h,L)=cntpart0(:,h,L)/partmass(L)
!calc dpos 
	    i=(h-1)*nav(KOM)
            do m=1,nav(KOM)
              dpos(:,i+m,kom)=pos(:,i+m,kom)-cntpart0(:,h,L)
!write(*,*) dpos(:,i+m,kom)
	    enddo
          enddo
	ENDDO
!stop

!	^^^ record periodic images of cntpart ^^^
        iser=0
 	do ix=-ncopy,+ncopy ; do iy=-ncopy,+ncopy ; do iz=-ncopy,+ncopy
        iser=iser+1
        do L=1,nparts
        KOM=kom_part(L)+1
 	do h=1,nmv(kom)
	  cntpart(:,h,L,iser)=cntpart0(:,h,L)  &
     &                        +ix*avec(:) + iy*bvec(:) + iz*cvec(:)
 	enddo ! h
	enddo; enddo; enddo
 	enddo ! kom

!	^^^ calc r between mass centers, and order parameter ^^^!
        DO L=1,nparts
        KOM=kom_part(L)+1
 	do hi=1,nmv(kom)
	do iser=1,isermax
 	do hj=1,nmv(kom)
 	    xij=cntpart(1,hj,L,iser)-cntpart(1,hi,L,isermin0)
 	    yij=cntpart(2,hj,L,iser)-cntpart(2,hi,L,isermin0)
   	    zij=cntpart(3,hj,L,iser)-cntpart(3,hi,L,isermin0)
	    rij2=xij**2+yij**2+zij**2
! cut-off
	  if(rij2==0d0) cycle
	  if(rij2.gt.Rmin2) cycle
            rij=sqrt(rij2)
! vector normal to planer backbone of mol 1.
            i=(hi-1)*nav(kom)
            ii1=i+ia(1) +1
            ij1=i+ja(1) +1
            veci1(:)=dpos(:,ij1,kom)-dpos(:,ii1,kom)
            ii2=i+ia(2) +1
            ij2=i+ja(2) +1
            veci2(:)=dpos(:,ij2,kom)-dpos(:,ii2,kom)
            call vecpro(veci1,veci2,veci3)
! vector normal to planer backbone of mol 2.
            j=(hj-1)*nav(kom)
            ji1=j+ia(1) +1
            jj1=j+ja(1) +1
            vecj1(:)=dpos(:,jj1,kom)-dpos(:,ji1,kom)
            ji2=j+ia(2) +1
            jj2=j+ja(2) +1
            vecj2(:)=dpos(:,jj2,kom)-dpos(:,ji2,kom)
            call vecpro(vecj1,vecj2,vecj3)
! debug
        costheta=0d0
        absvi=0d0; absvj=0d0
        do k=1,3
          costheta=costheta+veci3(k)*vecj3(k)
          absvi=absvi+veci3(k)**2
          absvj=absvj+vecj3(k)**2
        enddo
        costheta=costheta/sqrt(absvi)/sqrt(absvj)
        theta=acos(costheta)
        Sv=Sv+0.5d0*(3d0*costheta**2-1d0)
        itheta=int((theta/pi*180d0)/dtheta)+1
        hsttheta(itheta)=hsttheta(itheta)+1
!
	enddo ! hj
	enddo ! iser
	enddo ! hi
        enddo ! L
!debug
!       nsample=0
!       do i=1,ntheta
!         theta=(i-0.5d0)*dtheta-180d0
!         nsample=nsample+hsttheta(i)
!         write(*,*) theta,hsttheta(i)
!       enddo
!stop

        write(10,*) iflame, '-th flame analyzed/', nflame
        call flush(10)
      ENDDO ! iflame

      close(34)
      close(10)

!### average over flame ###!
!debug
       nsample=0
       do i=1,ntheta
         nsample=nsample+hsttheta(i)
!        write(*,*) theta,hsttheta(i)
       enddo
       do i=1,ntheta
         prbtheta(i)=hsttheta(i)/dble(nsample)
       enddo
       Sv=Sv/dble(nsample)

!### calc prob ###!
       dthetarad=dtheta/180d0*pi  ! rad
       prbs=0d0
       do i=1,ntheta
         theta=(i-0.5d0)*dtheta
         thetarad=abs(theta/180d0*pi) ! rad
         prbtheta(i)=prbtheta(i)/(2d0*pi*sin(thetarad)*dthetarad)
         prbs=prbs+prbtheta(i)
       enddo

!###  output  ###
       write(*,'(a,f10.5,i10)')'#', Sv, nsample
       do i=1,ntheta
         theta=(i-0.5d0)*dtheta
         prbtheta(i)=prbtheta(i)/prbs
         write(*,*) theta,prbtheta(i)
       enddo

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

!########################################################################
        subroutine vecpro(vecj1,vecj2,vecj3)
!########################################################################
        implicit none
        real(8) vecj1(3), vecj2(3), vecj3(3)

        vecj3(1)=vecj1(2)*vecj2(3) - vecj1(3)*vecj2(2)
        vecj3(2)=vecj1(3)*vecj2(1) - vecj1(1)*vecj2(3)
        vecj3(3)=vecj1(1)*vecj2(2) - vecj1(2)*vecj2(1)

        return
        end subroutine vecpro
!########################################################################
