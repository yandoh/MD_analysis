#ifort -O3 ../Rg_sa.F90 -o ../../exe/Rg
gfortran -O3 ../Rg_sa.F90 -o ../../exe/Rg
ln -sf npt120K1atm_mol.dcd ./DCD
../../exe/Rg > result.rg
