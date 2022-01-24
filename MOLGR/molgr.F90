!Copyright (c) 2022 Yoshimichi ANDOH
!Released under the MIT license
!https://opensource.org/licenses/mit-license.php

! === inputs ===
! (1) ./sys_info
! (2) ./3x3 
!
	implicit none
	integer(4),parameter::ncopy=1
	character(1)::cdum
	integer(4),allocatable::nav(:),nmv(:)
	integer(4)::h,i,j,k,l,m,n,komtot,kom,i0
	integer(4)::komi,komj,hi,hj
	integer(4)::iser,isermin1
 	integer(4),parameter :: isermax=(ncopy+1+ncopy)**3
	real(8)::xij,yij,zij,rij,rij2
	integer(4) :: totmol
	integer(4) :: totnp
	real(8) :: alpha,beta,gamma,box(3)
	real(8),allocatable :: cntmol0(:,:,:,:)
        integer(4) :: isermin0
	real(8),allocatable :: hst(:),gr(:)
	real(8) :: dr,rmax,rmax2
	integer(4) :: ir,ndr
	integer(4) :: molcount
	real(8) :: rho,radius,dvolume,volume,svolume
	real(8) :: pi
!dcd
	integer(4)::nflame,iflame,rflame
!dcd

	rmax=50d0
	rmax2=rmax*rmax
	dr=0.25d0
	ndr=int(rmax/dr)
	pi=dacos(-1d0)

	allocate(hst(ndr))
	allocate(gr(ndr))
	hst=0d0
	gr=0d0
	molcount=0
	volume=0d0

!### input sysinfo ###
	open(21,file='sys_info')
        read(21,*) ! skip
 	read(21,*) komtot
 	allocate( nav(komtot), nmv(komtot) )
        do kom=1,komtot
          read(21,*) nav(kom),nmv(kom)
          molcount=molcount+nmv(kom)
        enddo
	close(21)

	isermin0=(ncopy+1+ncopy)**3/2+1

        allocate(cntmol0(3,komtot,molcount,(ncopy+1+ncopy)**3))
	
!### input copied COM ###!
	read(*,*) cdum,nflame

      rflame=0
      DO iflame=1,nflame
        rflame=rflame+1

!	^^^ input cntmol ^^^
 	do kom=1,komtot
	do iser=1,isermax
 	do h=1,nmv(kom)
	  read(*,*) cntmol0(:,kom,h,iser)
 	enddo ! h
 	enddo ! kom
	enddo ! iser
	read(*,*) cdum,alpha,beta,gamma
	read(*,*) cdum,box(1:3)

        call cellvolume(box,alpha,beta,gamma,volume)
	svolume=svolume+volume

!	^^^ calc g(r) ^^^!
 	do komi=1,komtot
 	do hi=1,nmv(komi)
	do iser=1,isermax
 	do komj=1,komtot
 	do hj=1,nmv(komj)
 	    xij=cntmol0(1,komj,hj,iser)-cntmol0(1,komi,hi,isermin0)
 	    yij=cntmol0(2,komj,hj,iser)-cntmol0(2,komi,hi,isermin0)
   	    zij=cntmol0(3,komj,hj,iser)-cntmol0(3,komi,hi,isermin0)
	    rij2=xij**2+yij**2+zij**2
	  if(rij2==0d0) cycle
	  if(rij2.gt.rmax2) cycle
	    rij=sqrt(rij2)
	    ir=dint(rij/dr)+1
	    hst(ir)=hst(ir)+1d0
	enddo ! h
	enddo ! kom
	enddo ! iser
	enddo ! h
	enddo ! kom

      ENDDO ! iflame


!### average over flame ###!
      hst=hst/dble(rflame)
      svolume=svolume/dble(rflame)

!### calc gr ###!
      rho=molcount/svolume !(HH(1,1)*HH(2,2))

      do ir=1,ndr
 	radius=dr*(ir-0.5d0)
	dvolume=(4d0*pi*radius**2)*dr
 	gr(ir)=hst(ir)/molcount/(dvolume*rho)  
      enddo

!###  output  ###
      write(*,'(a1,2e22.15,i5,i7)') '#', rho, svolume, ndr, rflame
      do ir=1,ndr
 	radius=dr*(ir-0.5d0)
	write(*,'(f12.5,4e20.10)') radius,gr(ir)
      enddo

      stop
      end

!############################################################
        subroutine cellvolume(box,alpha,beta,gamma,volume)
        implicit none
        real(8)::box(3),alpha,beta,gamma,volume
	real(8),parameter::pi=dcos(-1d0)
	real(8)::av(3),bv(3),cv(3)
        real(8)::bvcv(3), avdotbvcv
        integer(4)::k

        call createvectors(box,alpha,beta,gamma,av,bv,cv)

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

        return
        end
