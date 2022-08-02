ifort -O3 ../atomgr_sa.F90 -o ../atomgr
#gfortran -O3 ../atomgr_sa.F90 -o ../atomgr
ln -sf water_npt_10ps.dcd ./DCD
../atomgr > result.gr
