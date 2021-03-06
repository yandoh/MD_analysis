!Copyright (c) 2022 Yoshimichi ANDOH
!Released under the MIT license
!https://opensource.org/licenses/mit-license.php

!!
!! Files to be required 
!! (1) ./sys_info
!! (2) ./massinfo.mdff  
!! (3) *.mdmntr files  (stdin)
!!
        implicit none
!! Number of data lines in .mdmntr file
        integer(4),parameter::nline=500000        
!! The Avogadro constant
        real(8),parameter :: avonum=6.0221409d+23 
!! Unit conversion
        real(8),parameter :: atm=1.013d+5
!
        integer(4),allocatable :: nav(:), nmv(:)
        real(8),allocatable :: molmass(:)
        real(8),allocatable :: vmole(:), density(:), kjmol(:)
!
        integer(4) :: kom, komtot
        integer(4) :: molcount
        real(8) :: massin
        integer(4)::nfiles
        integer(4)::h,i,j,nlall
        real(8) :: datin(20)
        real(8) :: avedat(20),avedat2(20)
        real(8) :: dev(20), std(20)

!### input sysinfo ###
!! 系の分子構成ファイルの読み込み
        molcount=0
	open(21,file='sys_info',status='old')
        read(21,*) ! skip
 	read(21,*) komtot
 	allocate( nav(komtot), nmv(komtot) )
        do kom=1,komtot
          read(21,*) nav(kom),nmv(kom)
          molcount=molcount+nmv(kom)
        enddo
	close(21)

	allocate( vmole(komtot), density(komtot), kjmol(komtot) )
        allocate(molmass(komtot))

        do kom=1,komtot
          vmole(kom)=dble(nmv(kom))/avonum
          kjmol(kom)=1d0/vmole(kom)*1d-3
	enddo
        
!### input massinfo ###
!! 質量情報ファイルの読み込み
        open(1,file='./massinfo.mdff',status='old')  
        molmass=0d0
        do kom=1,komtot
	read(1,*) ! skip <mass>
        do h=1,nav(kom)
          read(1,*) massin
          molmass(kom)=molmass(kom)+massin
        enddo
	read(1,*) ! skip </mass>
        enddo
        close(1)

        avedat=0d0
        avedat2=0d0
        nlall=0
        nfiles=0
1       continue
        read(*,*,end=99)  !! .mdmntr のヘッダーを読み飛ばす
        read(*,*)
        read(*,*)
        read(*,*)
        read(*,*)
        read(*,*)
        do i=1,nline
            read(*,*,end=99) j, datin(1:20)
!>>unit conversion
            datin(7) =datin(7) *1d+30
            datin(9) =datin(9) *1d+10
            datin(10)=datin(10)*1d+10
            datin(11)=datin(11)*1d+10
            datin(8)=datin(8)/atm
            datin(15:20)=datin(15:20)/atm
!<<unit conversion
            avedat(:)=avedat(:)+datin(:)
            avedat2(:)=avedat2(:)+datin(:)**2
            nlall=nlall+1
        enddo
        nfiles=nfiles+1
        goto 1            !! 次のファイルへ移動
99      continue

!average
        avedat=avedat/dble(nlall)
        avedat2=avedat2/dble(nlall)
        dev=avedat2-avedat**2
        std=sqrt(dev)

        do kom=1,komtot
          density(kom)=vmole(kom)*molmass(kom)/avedat(7) * 1d-6 * 1d+30
        enddo

        write(*,2) 'Tlines=',nlall
        write(*,2) 'NFiles=',nfiles
        do kom=1,komtot
          write(*,4) 'Density=',density(kom), ' [g/cm3], kom=', kom
        enddo
        write(*,3) 'Volume =',avedat(7), ' [A3]    ',  std(7)
!
        write(*,3) '|a|    =',avedat(9), ' [A]     ',   std(9)
        write(*,3) '|b|    =',avedat(10),' [A]     ',   std(10) 
        write(*,3) '|c|    =',avedat(11),' [A]     ',   std(11) 
        write(*,3) 'alpha  =',avedat(12),' [deg]   ', std(12)
        write(*,3) 'beta   =',avedat(13),' [deg]   ', std(13)
        write(*,3) 'gamma  =',avedat(14),' [deg]   ', std(14)
!
        write(*,3) 'Temp   =',avedat(6), ' [K]     ',   std(6)
        write(*,3) 'Press  =',avedat(8), ' [atm]   ', std(8)
        write(*,3) 'Pxx    =',avedat(15),' [atm]   ', std(15)
        write(*,3) 'Pyy    =',avedat(16),' [atm]   ', std(16)
        write(*,3) 'Pzz    =',avedat(17),' [atm]   ', std(17)
        write(*,3) 'Pxy    =',avedat(18),' [atm]   ', std(18)
        write(*,3) 'Pxz    =',avedat(19),' [atm]   ', std(19)
        write(*,3) 'Pyz    =',avedat(20),' [atm]   ', std(20)
        write(*,3) 'Pot_E  =',avedat(2)*kjmol,' [kJ/mol]',std(2)*kjmol
2       format(a,i23,a)
!3       format(a,es23.15,a,es23.15)
3       format(a,f23.6,a,f23.6)
4       format(a,f23.6,a,i5)

        stop
        end
