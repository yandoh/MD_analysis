        implicit none
        real(8),parameter :: atm=1.013d+5
        integer(4),parameter::nline=500000        !! Number of data lines in .mdmntr file
        real(8),parameter :: avonum=6.0221409d+23 !! The Avogadro constant
!
        integer(4),allocatable :: nav(:), nmv(:)
        real(8),allocatable :: molmass(:)
        real(8),allocatable :: vmole(:), density(:), jmol(:)
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
        molcount=0
	open(21,file='sys_info')
        read(21,*) ! skip
 	read(21,*) komtot
 	allocate( nav(komtot), nmv(komtot) )
        do kom=1,komtot
          read(21,*) nav(kom),nmv(kom)
          molcount=molcount+nmv(kom)
        enddo
	close(21)

	allocate( vmole(komtot), density(komtot), jmol(komtot) )
        allocate(molmass(komtot))

        do kom=1,komtot
          vmole(kom)=dble(nmv(kom))/avonum
          jmol(kom)=1d0/vmole(kom)*1d-3
	enddo
        
!### input massinfo ###
        open(1,file='./massinfo.mdff')    !! 質量情報ファイルの読み込み
        molmass=0d0
        do kom=1,komtot
	read(1,*)
        do h=1,nav(kom)
          read(1,*) massin
          molmass(kom)=molmass(kom)+massin
        enddo
	read(1,*)
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
            datin(8)=datin(8)/atm
            datin(15:20)=datin(15:20)/atm
!<<unit conversion
            avedat(:)=avedat(:)+datin(:)
            avedat2(:)=avedat2(:)+datin(:)**2
!               write(*,*) i,j
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
          density(kom)=vmole(kom)*molmass(kom)/avedat(7) * 1d-6
        enddo

        write(*,2) 'Tlines=',nlall
        write(*,2) 'NFiles=',nfiles
        do kom=1,komtot
          write(*,4) 'Densit=',density(kom), ' [g/cm3], kom=', kom
        enddo
        write(*,3) 'Volume=',avedat(7), std(7),   ' [m3]'
!
        write(*,3) '|a|   =',avedat(9), std(9),   ' [m]'
        write(*,3) '|b|   =',avedat(10), std(10), ' [m]'
        write(*,3) '|c|   =',avedat(11), std(11), ' [m]'
        write(*,3) 'alpha =',avedat(12), std(12), ' [deg]'
        write(*,3) 'beta  =',avedat(13), std(13), ' [deg]'
        write(*,3) 'gamma =',avedat(14), std(14), ' [deg]'
!
        write(*,3) 'Temper=',avedat(6), std(6),   ' [K]'
        write(*,3) 'Press =',avedat(8), std(8),   ' [atm]'
        write(*,3) 'Pxx   =',avedat(15), std(15), ' [atm]' 
        write(*,3) 'Pyy   =',avedat(16), std(16), ' [atm]'
        write(*,3) 'Pzz   =',avedat(17), std(17), ' [atm]'
        write(*,3) 'Pxy   =',avedat(18), std(18), ' [atm]'
        write(*,3) 'Pxz   =',avedat(19), std(19), ' [atm]'
        write(*,3) 'Pyz   =',avedat(20), std(20), ' [atm]'
        write(*,3) 'Pot_E =',avedat(2)*jmol, std(2)*jmol, ' [kJ/mol]'
2       format(a,i23,a)
3       format(a,2es23.15,a)
4       format(a,f23.10,a,i5)

        stop
        end
