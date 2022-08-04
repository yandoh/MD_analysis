COMPILER=gfortran
$COMPILER -O3 ../atomgr_sa.F90 -o ../../exe/atomgr
ln -sf water_npt_10ps.dcd ./DCD
../../exe/atomgr > result.gr
