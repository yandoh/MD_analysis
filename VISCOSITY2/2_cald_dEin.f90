!Copyright (c) 2022 Yoshimichi ANDOH
!Released under the MIT license
!https://opensource.org/licenses/mit-license.php

        implicit none
        real(8),parameter :: dt=2d-15
        integer(4) :: nl, nline, nomp, i
        real(8) :: volume, vkbt, temp
        character(1)::dum
        real(8),allocatable::time(:)
        real(8),allocatable::val(:,:)
        real(8),allocatable::nondiag(:)
!
        real(8)::diff(4:6), diffnond
        real(8)::dval(4:6), dnondiag

!###input
        read(*,*) dum, volume, vkbt, nline, temp, nomp
        write(1,*) nline
        allocate(time(nline),val(6,nline)) 
        allocate(nondiag(nline))
	nondiag=0d0
!stop
        do i=1,nline
          read(*,*) time(i), val(1:6,i)
          nondiag(i)=sum(val(4:6,i))/3d0
        enddo
!stop
        nl=0
        dval=0d0
        dnondiag=0d0
        write(*,1) dum, volume, vkbt, nline, temp, nomp
!###take difference
        do i=1,nline-1
          nl=nl+1
!
          diff(4:6)=(val(4:6,i+1)-val(4:6,i))/dt
	  diffnond=(nondiag(i+1)-nondiag(i))/dt
!
          dval(4:6)=dval(4:6)+diff(4:6)
          dnondiag=dnondiag+diffnond
!
          write(*,2) time(i), dnondiag/dble(nl), dval(4:6)/dble(nl)
        enddo
1       format(a1,2es23.15,i10,f10.2,i5)
2       format(f15.6,8es23.15)

        stop
        end
