!Copyright (c) 2022 Yoshimichi ANDOH
!Released under the MIT license
!https://opensource.org/licenses/mit-license.php

!!
!! Files to be required 
!! (1) ./sys_info
!! (2) ./seginfo.mdff  
!! (3) ./massinfo.mdff  
!! (4) *.dcd file 
!!
	implicit none
	real(8),allocatable::pos(:,:), xyz(:,:),vel(:,:)
	character(4),allocatable::header(:)
	character(4),allocatable::atmname(:)
	character(4),allocatable::resname(:)
	integer(4),allocatable::resnum(:)
	integer(4),allocatable::nav(:),nmv(:),nsegv(:)
	integer(4)::k0,L0,j0
	integer(4)::h,h0,i,k,j,l,m,n,komtot,kom,navmax
	real(8)::box(3), wa(3)
	real(8)::alpha,beta,gamma
        integer(4) :: numseg, idum, idum2, idum3
	integer(4),allocatable :: memseg(:),memsegin(:)
	integer(4),allocatable :: sumseg(:),sumsegin(:)
	integer(4),allocatable :: topseg(:),topsegin(:)
	real(8),allocatable :: cntseg(:,:), scntseg(:,:)
	real(8),allocatable :: dseg(:,:)
	real(8),allocatable :: cntmol(:,:), dmol(:,:)
	real(8),allocatable :: scntmol(:,:)
	real(8),allocatable :: mass(:),massin(:)
	real(8) :: smass
	real(8),allocatable :: molmass(:)
	real(8) :: sum(3),ave(3)
	real(8) :: pai, radi
	real(8) :: av(3),bv(3),cv(3), HH(3,3), rH(3,3)
	real(8) :: sxij,syij,szij
	integer(4) :: isum,isumplus
	integer(4) :: ksum
	integer(4) :: atmsum, molnum
	integer(4) :: imol,jmol
	integer(4) :: hsum,totmol,totatm
	integer(4) :: imem
	integer(4),allocatable :: seg2mol(:), seg2kom(:)
!dcd
      character(len=4)::Aname
      integer(4) :: nstr, nptot, iflame
      integer(4),dimension(20)::icntrl
      real(8),dimension(6)::cellstr
      real(4),allocatable::flpx(:),flpy(:),flpz(:)
      real(8)::pi
      Aname='CORD'
      icntrl=0
      nstr=0
      pi=acos(-1d0)
!dcd

!### input sysinfo ###
  	open(21,file='sys_info',status='old')
        read(21,*)
 	read(21,*) komtot
  	allocate( nav(komtot), nmv(komtot), nsegv(komtot) )
        do kom=1,komtot
 	  read(21,*) nav(kom),nmv(kom)
        enddo
        close(21)

	totmol=0
	totatm=0
	navmax=0
	do kom=1,komtot
	  totmol=totmol+nmv(kom)
	  totatm=totatm+nmv(kom)*nav(kom)
	  navmax=max(nav(kom),navmax)
	enddo
	allocate( cntmol(3,totmol), dmol(3,totatm) )
	allocate( scntmol(3,totmol) )
	allocate( molmass(komtot) )
 	allocate( pos(3,totatm),xyz(3,totatm) )
 	allocate( vel(3,totatm) )
	n=totatm
 	allocate(header(n))
 	allocate(atmname(n))
 	allocate(resname(n))
 	allocate(resnum(n))
 	allocate( mass(totatm) )
 	allocate( massin(navmax) )

!       write(1,*) totatm
!       stop

!### input seginfo.mdff ###
	open(20,file='seginfo.mdff',status='old')
	numseg=0
	do kom=1,komtot
	  read(20,*) ! skip tag <segments>
	  read(20,*) idum   !; write(*,*) 'idum=',idum
	  nsegv(kom)=idum
	  numseg=numseg+idum*nmv(kom)
	  do i=1,idum
	    read(20,*) ! skip tag <segment>
	    read(20,*) ! ID
	    read(20,*) idum2 !; write(*,*) 'idum2=',idum2
	    read(20,*) ! skip tag <atom>
	    do j=1,idum2
	      read(20,*) ! skip num
	    enddo
	    read(20,*) ! skip tag </atom>
	    read(20,*) ! skip tag </segment>
	  enddo
	  read(20,*) ! skip tag </segments>
	enddo ! kom
	close(20)

 	allocate(memseg(totatm),sumseg(numseg),topseg(numseg))
 	allocate(memsegin(navmax),sumsegin(2000),topsegin(2000))
 	allocate(cntseg(3,numseg),scntseg(3,numseg))
 	allocate(dseg(3,totatm))
 	allocate(seg2mol(numseg))
 	allocate(seg2kom(numseg))
 	scntseg=0d0


	open(20,file='seginfo.mdff',status='old')

	k=0;L=0
	DO kom=1,komtot

	  L0=0
	  !== input data ===!
	  read(20,*) ! skip tag <segments>
	  read(20,*) idum ! nsegment in kom-th mole
	  do k0=1,idum ! local segment number
	    read(20,*) ! skip tag <segment>
	    read(20,*) ! skip ID
	    read(20,*) idum2 !; write(*,*) idum2
	    sumsegin(k0)=idum2
	    read(20,*) ! skip tag <atom>
	    do j0=1,idum2
	      L0=L0+1  ! local atom number
	      read(20,*) idum3 !; write(*,*) '   ',idum3+1
            memsegin(L0)=idum3+1
	      if(j0.eq.1) topsegin(k0)=idum3+1
	    enddo
	    read(20,*) ! skip tag </atom>
	    read(20,*) ! skip tag </segment>
	  enddo
	  read(20,*) ! skip tag </segments>

	  !== copy data ===!
	  isum=0;ksum=0
	  if(kom.ge.2)then
	    do j=1,kom-1
	      isum=isum+nmv(j)*nav(j)
	      ksum=ksum+nmv(j)*nsegv(j)
	    enddo ! j
	  endif

	  do h=1,nmv(kom)
	  L0=0
	  do k0=1,nsegv(kom)
	    sumseg(k0+(h-1)*nsegv(kom)+ksum)=sumsegin(k0)
	    do j0=1,sumsegin(k0)
	      L0=L0+1  ! local atom number
	      memseg(L0+(h-1)*nav(kom)+isum)=memsegin(L0)+(h-1)*nav(kom)+isum
	      if(j0==1) topseg(k0+(h-1)*nsegv(kom)+ksum)=memseg(L0+(h-1)*nav(kom)+isum)
	    enddo ! j0
	  enddo ! k0
	  enddo ! h

	ENDDO ! kom
	close(20)


      !--> create seg2kom, seg2mol <--
	isum=0
	isumplus=0
	kom=1
	atmsum=0
	molnum=1
	do i=1,numseg
	
	  atmsum=atmsum+sumseg(i)
	  if( atmsum.le.nav(kom) )then
	    seg2mol(i)=molnum
	    if( atmsum.eq.nav(kom) )then
	       molnum=molnum+1
	       atmsum=0
	    endif
	  endif
	
	  isum=isum+sumseg(i) !+isumplus
	  if(isum.le.nav(kom)*nmv(kom)+isumplus)then
	    seg2kom(i)=kom
	    if(isum.eq.nav(kom)*nmv(kom)+isumplus)then
	       isumplus=isumplus+nav(kom)*nmv(kom)
	       kom=kom+1
	    endif
	  endif

	enddo

!### input massinfo ###
 	open(30,file='massinfo.mdff',status='old')  !! copied from mdff
	do kom=1,komtot
	  ! read !
	  read(30,*)
	  do m=1,nav(kom)
	    read(30,*) massin(m) !, cdum, atmname(m)
            call detect_atomkind(massin(m),atmname(m))
!           write(*,*) atmname(m)
	  enddo ! m
	  read(30,*)
	  ! set value !
	  isum=0
	  if(kom.ge.2)then
	    do j=1,kom-1
	      isum=isum+nmv(j)*nav(j)
	    enddo ! j
	  endif
	  do h=1,nmv(kom)
	  i=(h-1)*nav(kom)+isum
	  do m=1,nav(kom)
	    mass(i+m)=massin(m)
	  enddo ! m
	  enddo ! h
	enddo ! kom
 	close(30)

!### input .dcd file ###
	open(34,file='./DCD',form='unformatted',status='old')
	read(34) Aname,icntrl
	read(34) nstr
	read(34) nptot
        allocate( flpx(nptot),flpy(nptot),flpz(nptot) )

!### output .dcd file ###
	open(35,file='./new.dcd',form='unformatted',status='replace')
	write(35) Aname,icntrl
	write(35) nstr
	write(35) nptot

  	write(*,*) nptot,icntrl(1)

	DO iflame=1,icntrl(1)
	if(mod(iflame,100)==0) write(*,*) 'Converting ',iflame,'-th flame'
	read(34) cellstr !*1e+10
	read(34) (flpx(i),i=1,nptot)
	read(34) (flpy(i),i=1,nptot)
	read(34) (flpz(i),i=1,nptot)

	write(35) cellstr !*1e+10

        do i=1,nptot
          xyz(1,i)=flpx(i)
          xyz(2,i)=flpy(i)
          xyz(3,i)=flpz(i)
        enddo
 
      box(1)=cellstr(1)   ! |a|
      gamma =cellstr(2)   ! degree
      box(2)=cellstr(3)   ! |b|
      beta  =cellstr(4)   ! degree
      alpha =cellstr(5)   ! degree
      box(3)=cellstr(6)   ! |c|

	!write(*,*) box
	!write(*,*) alpha,beta,gamma

!### geometric center of segment ###
	cntseg=0d0
	L=0
	do i=1,numseg
	  wa=0d0
	  smass=0d0
	  do j=1,sumseg(i)
	    L=L+1
	    wa(:)=wa(:)+xyz(:,memseg(L))*mass(L)
	    smass=smass+mass(L)
	  enddo ! j
	  cntseg(:,i)=wa(:)/smass
	enddo ! i

	L=0
	do i=1,numseg
	do j=1,sumseg(i)
	L=L+1
	do k=1,3
	  dseg(k,memseg(L))=  &
              xyz(k,memseg(L))-cntseg(k,i)  ! r - r_{segcnt}
  	  if( abs(dseg(k,memseg(L))).ge.min(box(1),box(2),box(3))/2d0 ) stop 'ERROR'
	enddo
	enddo
	enddo

!### apply segment based PBC ###
        pai=dacos(-1d0)
        radi=pai/180d0
        alpha=alpha*radi
        beta =beta *radi
        gamma=gamma*radi

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
        HH(2,1)=av(2)
        HH(3,1)=av(3)
        HH(1,2)=bv(1)
        HH(2,2)=bv(2)
        HH(3,2)=bv(3)
        HH(1,3)=cv(1)
        HH(2,3)=cv(2)
        HH(3,3)=cv(3)
        call reverse(HH,rH)

!	^^^ cntseg 2 scntseg ^^^
	do i=1,numseg
	do k=1,3
	  sum(k)=0d0
	  do j=1,3
	    sum(k)=sum(k)+rH(k,j)*cntseg(j,i)
	  enddo
	  scntseg(k,i)=sum(k)
	
 	  if (scntseg(k,i).gt.+0.5d0) then
 	! do while (scntseg(k,i).gt.+0.5d0)
 	    scntseg(k,i)=scntseg(k,i)-1d0
        ! enddo
 	! do while (scntseg(k,i).le.-0.5d0)
 	  elseif (scntseg(k,i).le.-0.5d0) then
 	    scntseg(k,i)=scntseg(k,i)+1d0
 	! enddo
 	  endif
	
	enddo
	enddo

!### recovery of molecule shape ###
	sum=0d0
	imem=0
	ave=0d0
	
	do i=1,numseg-1
	imol=seg2mol(i  )
	jmol=seg2mol(i+1)
	
	if(jmol.eq.imol)then
	
	  do k=1,3
	    sum(k)=sum(k)+scntseg(k,i  )
	  enddo
	  imem=imem+1
	  ave=sum/dble(imem)
	
	  sxij=scntseg(1,i+1)-scntseg(1,i  )
	! sxij=scntseg(1,i+1)-ave(1)
	  do while (sxij.gt.+0.5d0)
	    scntseg(1,i+1)=scntseg(1,i+1)-1d0
	!   sxij=scntseg(1,i+1)-ave(1)
	    sxij=scntseg(1,i+1)-scntseg(1,i  )
	  enddo
	! elseif(sxij.le.-0.5d0)then
	  do while (sxij.le.-0.5d0)
	    scntseg(1,i+1)=scntseg(1,i+1)+1d0
	!   sxij=scntseg(1,i+1)-ave(1)
	    sxij=scntseg(1,i+1)-scntseg(1,i  )
	  enddo
	
	  syij=scntseg(2,i+1)-scntseg(2,i  )
	! syij=scntseg(2,i+1)-ave(2)
	! if(syij.gt.+0.5d0)then
	  do while (syij.gt.+0.5d0)
	    scntseg(2,i+1)=scntseg(2,i+1)-1d0
	!   syij=scntseg(2,i+1)-ave(2)
	    syij=scntseg(2,i+1)-scntseg(2,i  )
	  enddo
	! elseif(syij.le.-0.5d0)then
	  do while (syij.le.-0.5d0)
	    scntseg(2,i+1)=scntseg(2,i+1)+1d0
	!   syij=scntseg(2,i+1)-ave(2)
	    syij=scntseg(2,i+1)-scntseg(2,i  )
	  enddo
	
	  szij=scntseg(3,i+1)-scntseg(3,i  )
	! szij=scntseg(3,i+1)-ave(3)
	  do while (szij.gt.+0.5d0)
	! if(szij.gt.+0.5d0)then
	    scntseg(3,i+1)=scntseg(3,i+1)-1d0
	!   szij=scntseg(3,i+1)-ave(3)
	    szij=scntseg(3,i+1)-scntseg(3,i  )
	  enddo
	  do while (szij.le.-0.5d0)
	! elseif(szij.le.-0.5d0)then
	    scntseg(3,i+1)=scntseg(3,i+1)+1d0
	!   szij=scntseg(3,i+1)-ave(3)
	    szij=scntseg(3,i+1)-scntseg(3,i  )
	  enddo
	
	else
	
	  do k=1,3
	    sum(k)=scntseg(k,i+1)
	  enddo
	  imem=1  ! reset imem
	  ave=sum/dble(imem)
	
	endif
	enddo


!### recovery of real cntseg ###
	do i=1,numseg
	do k=1,3
	  sum(k)=0d0
	  do j=1,3
	    sum(k)=sum(k)+HH(k,j)*scntseg(j,i)
	  enddo
	  cntseg(k,i)=sum(k)
	enddo
	enddo

!### new-xyz ###
	L=0
	do i=1,numseg
	do j=1,sumseg(i)
	L=L+1
	do k=1,3
	  xyz(k,memseg(L))=cntseg(k,i)+dseg(k,memseg(L))
	enddo
	enddo
	enddo

!! molecule shape was recovered !!



!### calc. COM of molecule ###
	DO KOM=1,komtot
	  
	  isum=0
	  hsum=0
	  if(KOM.ge.2)then
	    do L=1,KOM-1
	      isum=isum+nmv(L)*nav(L)
	      hsum=hsum+nmv(L)
	    enddo
	  endif
	
	smass=0d0
!     --- molecule mass ---
      do m=1,nav(KOM) 
	  smass=smass+mass(isum+m)
      enddo
	molmass(KOM)=smass

!     --- COM of each molecule ---
	DO h=1,nmv(KOM)
	  i=(h-1)*nav(KOM)+isum
	  sum(:)=0d0
	  do k=1,3
	  do m=1,nav(KOM)
	    sum(k)=sum(k)+xyz(k,i+m)*mass(i+m)
	  enddo
	  enddo  ! k
	  sum(:)=sum(:)/molmass(KOM)
	  cntmol(:,h+hsum)=sum(:)
	  do k=1,3
	  do m=1,nav(KOM)
	    dmol(k,i+m)=xyz(k,i+m)-cntmol(k,h+hsum)
	  enddo
	  enddo  ! k
	ENDDO  ! h
	

	ENDDO

      call  PBC_cntmol(komtot,totmol,cntmol,rH,HH,nmv)
!     call  reset_z0(komtot,totmol,cntmol,rH,HH,molmass,nmv,nsolv)

!=== recovery of pos ===!
        DO KOM=1,komtot
        
          isum=0
          hsum=0
          if(KOM.ge.2)then
            do L=1,KOM-1
              isum=isum+nmv(L)*nav(L)
              hsum=hsum+nmv(L)
            enddo
          endif

!       === pos ===
	  do h=1,nmv(KOM)
	  i=(h-1)*nav(KOM)+isum
	  do m=1,nav(KOM)
	  do k=1,3
	    pos(k,i+m)=cntmol(k,h+hsum)+dmol(k,i+m)
	  enddo ! k
	  enddo ! m
	  enddo ! h

	ENDDO ! KOM

!99	continue
!### output pdb with TER and END ###
!       do kom=1,komtot
!	  isum=0
!	  if(kom.ge.2)then
!	    do j=1,kom-1
!	      isum=isum+nmv(j)*nav(j)
!	    enddo
!	  endif

!       do h=1,nmv(kom)
!       i=(h-1)*nav(kom)+isum
!       do m=1,nav(kom)
!         write(*,111) header(i+m), i+m, &
!         atmname(i+m), resname(i+m), resnum(i+m), (pos(k,i+m),k=1,3)
!       enddo
!       write(*,'(a3)') 'TER'
!       enddo
!       enddo
!	write(*,'(a3)') 'END'
!	write(*,*) alpha,beta,gamma
!	write(*,*) box
!	write(*,*) n
!111	format(a4,i7,1x,a4,1x,a4,i5,4x,3f8.3)

        do i=1,nptot
          flpx(i)=pos(1,i)
          flpy(i)=pos(2,i)
          flpz(i)=pos(3,i)
        enddo
!### output .dcd ###
	write(35) (flpx(i),i=1,nptot)
	write(35) (flpy(i),i=1,nptot)
	write(35) (flpz(i),i=1,nptot)

99	continue

      if(iflame==icntrl(1))then
!### output last.pdb ###
 	open(31,file='last.pdb',status='replace')
        write(31,'(a6,3f9.3,3f7.3)') 'CRYST1', &
        box(1),box(2),box(3),alpha/radi,beta/radi,gamma/radi
	h0=0
 	do kom=1,komtot
 	  isum=0
 	  if(kom.ge.2)then
 	    do j=1,kom-1
 	      isum=isum+nmv(j)*nav(j)
 	    enddo
 	  endif
 	  do h=1,nmv(kom)
          h0=h0+1
 	  i=(h-1)*nav(kom)+isum
 	  do m=1,nav(kom)
 	    write(31,120) "ATOM", i+m, &
                atmname(m), "    ", h0, flpx(i+m),flpy(i+m),flpz(i+m) !, (pos(k,i+m),k=1,3)
 	  enddo
120	format(a4,i7,1x,a4,1x,a4,i5,4x,3f8.3)
 	write(31,'(a)') 'TER' ! skip TER
 	  enddo
 	enddo
 	close(31)
      endif
	
	ENDDO  ! iflame

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
!       det=HH(1,3)*( HH(2,1)*HH(3,2)-HH(2,2)*HH(3,1) )  &
!          +HH(2,3)*( HH(3,1)*HH(1,2)-HH(3,2)*HH(1,1) )  &
!          +HH(3,3)*( HH(1,1)*HH(2,2)-HH(1,2)*HH(2,1) )
        if(det.eq.0d0) stop 'det = 0'
!       write(*,*) det 
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

!       write(*,*)'HH='
!       do i=1,3
!       write(*,*) (HH(i,k),k=1,3)
!       enddo
!       write(*,*)'rH='
!       do i=1,3
!       write(*,*) (rH(i,k),k=1,3)
!       enddo

        return
        end

!########################################################################
        subroutine PBC_cntmol(komtot,totmol,cntmol,rH,HH,nmv)
!########################################################################
        implicit none
        integer(4)::KOM,hsum,L,h,k,j,komtot,totmol
        real(8)::cntmol(3,totmol)
        real(8)::scntmol(3,totmol)
        real(8)::rH(3,3),HH(3,3)
        integer(4)::nmv(komtot)
        real(8)::sum(3)

!       ^^^ svector of cntmol ^^^
	DO KOM=1,komtot
	
          hsum=0
          if(KOM.ge.2)then
            do L=1,KOM-1
              hsum=hsum+nmv(L)
            enddo
          endif

	  do h=1,nmv(KOM)
	  do k=1,3
	    sum(k)=0d0
	    do j=1,3
	      sum(k)=sum(k)+rH(k,j)*cntmol(j,h+hsum)
	    enddo ! j
	    scntmol(k,h+hsum)=sum(k)
	  enddo ! k
	  enddo ! h
	
!         ^^^ apply molecule-based PBC (all direction) ^^^
	  do h=1,nmv(KOM)
!
	! do while (scntmol(1,h+hsum).gt.+0.5d0)
	  if(scntmol(1,h+hsum).gt.+0.5d0)then
	    scntmol(1,h+hsum)=scntmol(1,h+hsum)-1d0
	! enddo
	! do while (scntmol(1,h+hsum).le.-0.5d0)
	  elseif(scntmol(1,h+hsum).le.-0.5d0)then
	    scntmol(1,h+hsum)=scntmol(1,h+hsum)+1d0
	! enddo
	  endif
!
	! do while (scntmol(2,h+hsum).gt.+0.5d0)
	  if(scntmol(2,h+hsum).gt.+0.5d0)then
	    scntmol(2,h+hsum)=scntmol(2,h+hsum)-1d0
	! enddo
	! do while (scntmol(2,h+hsum).le.-0.5d0)
	  elseif(scntmol(2,h+hsum).le.-0.5d0)then
	    scntmol(2,h+hsum)=scntmol(2,h+hsum)+1d0
	! enddo
	  endif
!
	! do while (scntmol(3,h+hsum).gt.+0.5d0)
	  if(scntmol(3,h+hsum).gt.+0.5d0)then
	    scntmol(3,h+hsum)=scntmol(3,h+hsum)-1d0
	! enddo
	! do while (scntmol(3,h+hsum).le.-0.5d0)
	  elseif(scntmol(3,h+hsum).le.-0.5d0)then
	    scntmol(3,h+hsum)=scntmol(3,h+hsum)+1d0
	! enddo
	  endif
!
	  enddo
	
	ENDDO

!       ^^^ recovery of cntmol & pos ^^^
        DO KOM=1,komtot
        
          hsum=0
          if(KOM.ge.2)then
            do L=1,KOM-1
              hsum=hsum+nmv(L)
            enddo
          endif
!         === cntmol ===
          do h=1,nmv(KOM)
          do k=1,3
            sum(k)=0d0
            do j=1,3
              sum(k)=sum(k)+HH(k,j)*scntmol(j,h+hsum)
            enddo ! j
            cntmol(k,h+hsum)=sum(k)
	  enddo ! k
	  enddo ! h

	ENDDO ! KOM

        return
        end

      subroutine  reset_z0(komtot,totmol,cntmol,rH,HH,molmass,nmv,nsolv)
        implicit none
        integer(4)::nsolv,komtot,totmol
        real(8)::cntmol(3,totmol)
        real(8)::molmass(totmol)
        integer(4)::nmv(komtot)
        real(8) :: HH(3,3), rH(3,3)
!
        integer(4)::KOM,hsum,L,h
        real(8)::blycnt(3)
        real(8)::totmss

!### calc. com of two bilayers ###!
        blycnt=0d0
        totmss=0d0
        DO KOM=1,komtot-nsolv
        
          hsum=0
          if(KOM.ge.2)then
            do L=1,KOM-1
              hsum=hsum+nmv(L)
            enddo
          endif

          do h=1,nmv(KOM)
        blycnt(:)=blycnt(:)+molmass(KOM)*cntmol(:,h+hsum)
        totmss=totmss+molmass(KOM)
          enddo !h

        ENDDO

        blycnt=blycnt/totmss

!### reset z0 ###!
        DO KOM=1,komtot

          hsum=0
          if(KOM.ge.2)then
            do L=1,KOM-1
              hsum=hsum+nmv(L)
            enddo
          endif
          do h=1,nmv(KOM)
        cntmol(3,h+hsum)=cntmol(3,h+hsum)-blycnt(3)
          enddo !h

        ENDDO

      call  PBC_cntmol(komtot,totmol,cntmol,rH,HH,nmv)

        return
      end

      subroutine detect_atomkind(massi,atmname)
      implicit none
      real(8)::massi
      character(4)::atmname

        if (dabs(massi-1.00800d0) .lt. 1.0d-2) then
           atmname = 'H   '
        else if (dabs(massi-12.01100d0) .lt. 1.0d-2) then
           atmname = 'C   '
        else if (dabs(massi-14.00700d0) .lt. 1.0d-2) then
           atmname = 'N   '
        else if (dabs(massi-15.99900d0) .lt. 1.0d-2) then
           atmname = 'O   '
        else if (dabs(massi-30.974000d0) .lt. 1.0d-2) then
           atmname = 'P   '
        else if (dabs(massi-32.06000d0) .lt. 1.0d-2) then
           atmname = 'S   '
        else if (dabs(massi-4.00260d0) .lt. 1.0d-2) then
           atmname = 'He  '
        else if (dabs(massi-20.17970d0) .lt. 1.0d-2) then
           atmname = 'Ne  '
!       else if (dabs(massi-40.08000d0) .lt. 1.0d-2) then
!          atmname = 9
!       else if (dabs(massi-65.37000d0) .lt. 1.0d-2) then
!          atmname = 10
!       else if (dabs(massi-55.84700d0) .lt. 1.0d-2) then
!          atmname = 11
        else if (dabs(massi-22.98977d0) .lt. 1.0d-2) then
           atmname = 'SOD '
        else if (dabs(massi-35.45000d0) .lt. 1.0d-2) then
           atmname = 'CLA '
        else if (dabs(massi-39.102000d0) .lt. 1.0d-2) then
           atmname = 'POT '
        else
           atmname = 'X   '
        endif
      
      end subroutine
