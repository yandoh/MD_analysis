#ifort -O3 ../SelfD_sa.F90 -o ../SelfD
gfortran -O3 ../SelfD_sa.F90 -o ../SelfD
ln -sf water_npt_10ps.dcd ./DCD
../SelfD > result.selfd
