!Copyright (c) 2022 Yoshimichi ANDOH
!Released under the MIT license
!https://opensource.org/licenses/mit-license.php

!
! === inputs ===
! (1) ./result.gr
!
	implicit none
	character(1) cdum
	real(8),allocatable::r(:),gr(:,:)
	real(8),allocatable::rho(:)
	integer(4) i
        integer(4) nparts, ndr, nflame
        real(8) volume, dr, pi
        real(8) rr,dv
        real(8),allocatable::dn(:),sn(:)
        real(8),allocatable::dKB(:),sKB(:)

	pi=acos(-1d0)

!### read g(r) ###
	read(*,*) cdum, volume, nparts, ndr, nflame
        allocate( rho(nparts) )
	read(*,*) cdum, rho(1:nparts)
	allocate(r(ndr),gr(ndr,nparts))
	allocate(dn(nparts),sn(nparts))
	allocate(dKB(nparts),sKB(nparts))
	do i=1,ndr
          read(*,*) r(i),gr(i,1:nparts)
	enddo
	dr=r(2)-r(1)

!### calc. coordinate number ###
        sn=0d0
        sKB=0d0
	do i=1,ndr
          rr=r(i)
          dv=4d0*pi*(rr**2)*dr
          dn(1:nparts)=rho(1:nparts)*dv*gr(i,1:nparts)
          sn(1:nparts)=sn(1:nparts)+dn(1:nparts)
          dKB(1:nparts)=(gr(i,1:nparts)-1d0)*dv
          sKB(1:nparts)=sKB(1:nparts)+dKB(1:nparts)
          write(*,'(f15.3,99f15.7)') rr, &
     &             gr(i,1:nparts), sn(1:nparts), sKB(1:nparts)
	enddo

	stop
	end
