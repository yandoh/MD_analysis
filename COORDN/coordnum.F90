!Copyright (c) 2022 Yoshimichi ANDOH
!Released under the MIT license
!https://opensource.org/licenses/mit-license.php

! Stand alone version of molgr.F90
!
! === inputs ===
! (1) ./result.gr
!
	implicit none
	character(1) cdum
	real(8),allocatable::r(:),gr(:,:)
	integer(4) i
        integer(4) nparts, ndr, nflame
        real(8) volume, rho, dr, pi
        real(8) rr,dv
        real(8),allocatable::dn(:),sn(:)

	pi=acos(-1d0)

!### read g(r) ###
	read(*,*) cdum, volume, nparts, ndr, nflame
	read(*,*) cdum, rho
	allocate(r(ndr),gr(ndr,nparts))
	allocate(dn(nparts),sn(nparts))
	do i=1,ndr
          read(*,*) r(i),gr(i,1:nparts)
	enddo
	dr=r(2)-r(1)

!### calc. coordinate number ###
        sn=0d0
	do i=1,ndr
          rr=r(i)
          dv=4d0*pi*(rr**2)*dr
          dn(1:nparts)=rho*dv*gr(i,1:nparts)
          sn(1:nparts)=sn(1:nparts)+dn(1:nparts)
          write(*,'(f15.3,99f15.7)') rr, &
     &             gr(i,1:nparts), sn(1:nparts)
	enddo

	stop
	end
