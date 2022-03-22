!Copyright (c) 2022 Yoshimichi ANDOH
!Released under the MIT license
!https://opensource.org/licenses/mit-license.php

        real(8), parameter :: dt=2e-15
!       integer(4), parameter :: imax=5000000  ! 10ns
!       integer(4), parameter :: imax=2500000  ! 5ns
        integer(4), parameter :: imax=1500000  ! 3ns
!       integer(4), parameter :: imax=1000000  ! 2ns
!       integer(4), parameter :: imax=500000   ! 1ns
        !
        integer(4) :: ndat
        integer(4) :: i,j,k
        character(1) :: dum
        real(8) :: ave(2)
        real(8) :: Ctin(6)
        real(8),allocatable :: Ct(:,:)
        integer(4),allocatable :: nhit(:)
        real(8) :: sekibun(6)
        real(8) :: sekibun_nond
        real(8) :: aved,avend
        real(8) :: maxaved,minaved
        real(8) :: maxavend,minavend

        sekibun=0
        sekibun_nond=0

!#  5.499431780322072E-26  1.358765895952437E-05   2500000
        read(*,*) dum, ave(1:2), ndat
        allocate(Ct(6,ndat),nhit(ndat))
        do i=1,ndat
          read(*,*) j, Ct(1:6,i), nhit(i)
!simple rule
         !sekibun(:)=sekibun(:)+dt*Ct(:,i)
         !sekibun_nond=sekibun_nond+dt*(Ct(4,i)+Ct(5,i)+Ct(6,i))
!trapezoidal rule
          if(i==1 .or. i==imax)then
            sekibun(:)=sekibun(:)+0.5d0*dt*Ct(:,i)
            sekibun_nond=sekibun_nond+0.5d0*dt*(Ct(4,i)+Ct(5,i)+Ct(6,i))
          else
            sekibun(:)=sekibun(:)+1.0d0*dt*Ct(:,i)
            sekibun_nond=sekibun_nond+1.0d0*dt*(Ct(4,i)+Ct(5,i)+Ct(6,i))
          endif
          write(*,'(f12.3,8es23.15)') j*dt*1e+12, &
      &        sekibun_nond/3d0,sekibun(4:6)
          if(i==imax) exit
        enddo

        sekibun_nond=sekibun_nond/3d0

        avend=(sekibun(4)+sekibun(5)+sekibun(6))/3d0
        maxavend=maxval(sekibun(4:6))
        minavend=minval(sekibun(4:6))
        write(*,4)'#imax=',imax
        write(*,3)'#>>  Viscosity(shear) / Pa*sec  <<<'
        write(*,2)'#=====  Ave.  =========', &
                  '#=====  Max.  =========', &
                  '#=====  Min.  =========|'
        write(*,1)'#',avend, maxavend, minavend

1       format(a,3es23.15)
2       format(3a)
3       format(a)
4       format(a,i10)

        stop
        end
