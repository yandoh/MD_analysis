#ifort -O3 ../partgr_sa.F90 -o ../../exe/partgr
gfortran -O3 ../partgr_sa.F90 -o ../../exe/partgr
ln -sf npt120K1atm_mol.dcd ./DCD
../../exe/partgr > result.rg
