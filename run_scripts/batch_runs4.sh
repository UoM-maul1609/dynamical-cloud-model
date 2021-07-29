#! /bin/bash

# batch runs - these are the runs using the sounding

# run 1: default
git checkout run_scripts/westbrook_illingworth_nml.in
git checkout run_scripts/westbrook_illingworth_pamm_nml.in 

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml.in 

mv /tmp/output_pc01.nc /tmp/output0001.nc


# run 2: no aerosol above
git checkout run_scripts/westbrook_illingworth_nml.in
git checkout run_scripts/westbrook_illingworth_pamm_nml.in 

sed -e "s/n_read(1,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(1,1:4)     = 250e6, 250e6, 250e6,0.1e6,/" run_scripts/westbrook_illingworth_pamm_nml.in > /tmp/namelist.tmp
sed -e "s/n_read(2,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(2,1:4)     = 250e6, 250e6, 250e6,0.1e6,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
sed -e "s/z_read(1:4)       = 0.,740,3260,4000/z_read(1:4)       = 0.,740,3700,3701/" /tmp/namelist.tmp2 > run_scripts/westbrook_illingworth_pamm_nml.in 

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml.in 

mv /tmp/output_pc01.nc /tmp/output0002.nc


# run 3: no aerosol above and mode 2 on
git checkout run_scripts/westbrook_illingworth_nml.in
git checkout run_scripts/westbrook_illingworth_pamm_nml.in 

sed -e "s/n_read(1,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(1,1:4)     = 250e6, 250e6, 250e6,0.1e6,/" run_scripts/westbrook_illingworth_pamm_nml.in > /tmp/namelist.tmp
sed -e "s/n_read(2,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(2,1:4)     = 250e6, 250e6, 250e6,0.1e6,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
sed -e "s/z_read(1:4)       = 0.,740,3260,4000/z_read(1:4)       = 0.,740,3700,3701/" /tmp/namelist.tmp2 > run_scripts/westbrook_illingworth_pamm_nml.in 

sed -e "s/nm1%mode2_ice_flag=0,/nm1%mode2_ice_flag=1,/" run_scripts/westbrook_illingworth_pamm_nml.in > /tmp/namelist.tmp
cp /tmp/namelist.tmp run_scripts/westbrook_illingworth_pamm_nml.in

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml.in 

mv /tmp/output_pc01.nc /tmp/output0003.nc


# run 4: no aerosol above and coll break-up on
git checkout run_scripts/westbrook_illingworth_nml.in
git checkout run_scripts/westbrook_illingworth_pamm_nml.in 

sed -e "s/n_read(1,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(1,1:4)     = 250e6, 250e6, 250e6,0.1e6,/" run_scripts/westbrook_illingworth_pamm_nml.in > /tmp/namelist.tmp
sed -e "s/n_read(2,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(2,1:4)     = 250e6, 250e6, 250e6,0.1e6,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
sed -e "s/z_read(1:4)       = 0.,740,3260,4000/z_read(1:4)       = 0.,740,3700,3701/" /tmp/namelist.tmp2 > run_scripts/westbrook_illingworth_pamm_nml.in 

sed -e "s/nm1%coll_breakup_flag1=0,/nm1%coll_breakup_flag1=1,/" run_scripts/westbrook_illingworth_pamm_nml.in > /tmp/namelist.tmp
cp /tmp/namelist.tmp run_scripts/westbrook_illingworth_pamm_nml.in

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml.in 

mv /tmp/output_pc01.nc /tmp/output0004.nc

# run 5: no aerosol above and lawson on
git checkout run_scripts/westbrook_illingworth_nml.in
git checkout run_scripts/westbrook_illingworth_pamm_nml.in 

sed -e "s/n_read(1,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(1,1:4)     = 250e6, 250e6, 250e6,0.1e6,/" run_scripts/westbrook_illingworth_pamm_nml.in > /tmp/namelist.tmp
sed -e "s/n_read(2,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(2,1:4)     = 250e6, 250e6, 250e6,0.1e6,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
sed -e "s/z_read(1:4)       = 0.,740,3260,4000/z_read(1:4)       = 0.,740,3700,3701/" /tmp/namelist.tmp2 > run_scripts/westbrook_illingworth_pamm_nml.in 

sed -e "s/nm1%lawson=.false./nm1%lawson=.true./" run_scripts/westbrook_illingworth_pamm_nml.in > /tmp/namelist.tmp
cp /tmp/namelist.tmp run_scripts/westbrook_illingworth_pamm_nml.in

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml.in 

mv /tmp/output_pc01.nc /tmp/output0005.nc


# run 6: no aerosol above and recycling off
git checkout run_scripts/westbrook_illingworth_nml.in
git checkout run_scripts/westbrook_illingworth_pamm_nml.in 

sed -e "s/n_read(1,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(1,1:4)     = 250e6, 250e6, 250e6,0.1e6,/" run_scripts/westbrook_illingworth_pamm_nml.in > /tmp/namelist.tmp
sed -e "s/n_read(2,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(2,1:4)     = 250e6, 250e6, 250e6,0.1e6,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
sed -e "s/z_read(1:4)       = 0.,740,3260,4000/z_read(1:4)       = 0.,740,3700,3701/" /tmp/namelist.tmp2 > run_scripts/westbrook_illingworth_pamm_nml.in 

sed -e "s/nm1%recycle=.true./nm1%recycle=.false./" run_scripts/westbrook_illingworth_pamm_nml.in > /tmp/namelist.tmp
cp /tmp/namelist.tmp run_scripts/westbrook_illingworth_pamm_nml.in

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml.in 

mv /tmp/output_pc01.nc /tmp/output0006.nc




# run 7: no aerosol above and stochastic nucleation 1e-8
git checkout run_scripts/westbrook_illingworth_nml.in
git checkout run_scripts/westbrook_illingworth_pamm_nml.in 

sed -e "s/n_read(1,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(1,1:4)     = 250e6, 250e6, 250e6,0.1e6,/" run_scripts/westbrook_illingworth_pamm_nml.in > /tmp/namelist.tmp
sed -e "s/n_read(2,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(2,1:4)     = 250e6, 250e6, 250e6,0.1e6,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
sed -e "s/z_read(1:4)       = 0.,740,3260,4000/z_read(1:4)       = 0.,740,3700,3701/" /tmp/namelist.tmp2 > run_scripts/westbrook_illingworth_pamm_nml.in 

sed -e "s/nm1%ice_nuc_flag=1,/nm1%ice_nuc_flag=2,/" run_scripts/westbrook_illingworth_pamm_nml.in > /tmp/namelist.tmp
sed -e "s/nm1%j_stochastic=0.5e-9,/nm1%j_stochastic=1e-8,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
cp /tmp/namelist.tmp2 run_scripts/westbrook_illingworth_pamm_nml.in

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml.in 

mv /tmp/output_pc01.nc /tmp/output0007.nc



# run 8: no aerosol above and stochastic nucleation 5e-8
git checkout run_scripts/westbrook_illingworth_nml.in
git checkout run_scripts/westbrook_illingworth_pamm_nml.in 

sed -e "s/n_read(1,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(1,1:4)     = 250e6, 250e6, 250e6,0.1e6,/" run_scripts/westbrook_illingworth_pamm_nml.in > /tmp/namelist.tmp
sed -e "s/n_read(2,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(2,1:4)     = 250e6, 250e6, 250e6,0.1e6,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
sed -e "s/z_read(1:4)       = 0.,740,3260,4000/z_read(1:4)       = 0.,740,3700,3701/" /tmp/namelist.tmp2 > run_scripts/westbrook_illingworth_pamm_nml.in 

sed -e "s/nm1%ice_nuc_flag=1,/nm1%ice_nuc_flag=2,/" run_scripts/westbrook_illingworth_pamm_nml.in > /tmp/namelist.tmp
sed -e "s/nm1%j_stochastic=0.5e-9,/nm1%j_stochastic=5e-8,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
cp /tmp/namelist.tmp2 run_scripts/westbrook_illingworth_pamm_nml.in

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml.in 

mv /tmp/output_pc01.nc /tmp/output0008.nc


# run 9: no aerosol above and stochastic nucleation 1e-7
git checkout run_scripts/westbrook_illingworth_nml.in
git checkout run_scripts/westbrook_illingworth_pamm_nml.in 

sed -e "s/n_read(1,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(1,1:4)     = 250e6, 250e6, 250e6,0.1e6,/" run_scripts/westbrook_illingworth_pamm_nml.in > /tmp/namelist.tmp
sed -e "s/n_read(2,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(2,1:4)     = 250e6, 250e6, 250e6,0.1e6,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
sed -e "s/z_read(1:4)       = 0.,740,3260,4000/z_read(1:4)       = 0.,740,3700,3701/" /tmp/namelist.tmp2 > run_scripts/westbrook_illingworth_pamm_nml.in 

sed -e "s/nm1%ice_nuc_flag=1,/nm1%ice_nuc_flag=2,/" run_scripts/westbrook_illingworth_pamm_nml.in > /tmp/namelist.tmp
sed -e "s/nm1%j_stochastic=0.5e-9,/nm1%j_stochastic=1e-7,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
cp /tmp/namelist.tmp2 run_scripts/westbrook_illingworth_pamm_nml.in

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml.in 

mv /tmp/output_pc01.nc /tmp/output0009.nc


# run 10: no aerosol above and mixing and advection of aerosol off
git checkout run_scripts/westbrook_illingworth_nml.in
git checkout run_scripts/westbrook_illingworth_pamm_nml.in 

sed -e "s/n_read(1,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(1,1:4)     = 250e6, 250e6, 250e6,0.1e6,/" run_scripts/westbrook_illingworth_pamm_nml.in > /tmp/namelist.tmp
sed -e "s/n_read(2,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(2,1:4)     = 250e6, 250e6, 250e6,0.1e6,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
sed -e "s/z_read(1:4)       = 0.,740,3260,4000/z_read(1:4)       = 0.,740,3700,3701/" /tmp/namelist.tmp2 > run_scripts/westbrook_illingworth_pamm_nml.in 

sed -e "s/!             if((n<2).or.(n>14)) then/             if((n<2).or.(n>14)) then/" sgm/subgrid_3d.f90 > /tmp/temp.f90
sed -e "s/!             endif/             endif/" /tmp/temp.f90 > /tmp/temp2.f90
cp /tmp/temp2.f90 sgm/subgrid_3d.f90

sed -e "s/!if((nqc<2).or.(nqc>(n_mode+1))) then/if((nqc<2).or.(nqc>(n_mode+1))) then/" driver_code.f90 > /tmp/temp.f90
sed -e "s/!endif/endif/" /tmp/temp.f90 > /tmp/temp2.f90
cp /tmp/temp2.f90 driver_code.f90

make NETCDFLIB=-L/usr/lib/x86_64-linux-gnu/ NETCDFMOD=/usr/include/ FFLAGS='-O3 -w -o' FFLAGS2='-w -O3 -o'


mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml.in 

mv /tmp/output_pc01.nc /tmp/output0010.nc


git checkout run_scripts/westbrook_illingworth_nml.in
git checkout run_scripts/westbrook_illingworth_pamm_nml.in 
git checkout sgm/subgrid_3d.f90
git checkout driver_code.f90


make NETCDFLIB=-L/usr/lib/x86_64-linux-gnu/ NETCDFMOD=/usr/include/ FFLAGS='-O3 -w -o' FFLAGS2='-w -O3 -o'




# run 11: no aerosol above and no radiation
git checkout run_scripts/westbrook_illingworth_nml.in
git checkout run_scripts/westbrook_illingworth_pamm_nml.in

sed -e "s/n_read(1,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(1,1:4)     = 250e6, 250e6, 250e6,0.1e6,/" run_scripts/westbrook_illingworth_pamm_nml.in > /tmp/namelist.tmp
sed -e "s/n_read(2,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(2,1:4)     = 250e6, 250e6, 250e6,0.1e6,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
sed -e "s/z_read(1:4)       = 0.,740,3260,4000/z_read(1:4)       = 0.,740,3700,3701/" /tmp/namelist.tmp2 > run_scripts/westbrook_illingworth_pamm_nml.in


sed -e "s/nm1%radiation=.true./nm1%radiation=.false./" run_scripts/westbrook_illingworth_pamm_nml.in > /tmp/namelist.tmp
cp /tmp/namelist.tmp run_scripts/westbrook_illingworth_pamm_nml.in

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml.in

mv /tmp/output_pc01.nc /tmp/output0011.nc


# run 12: no aerosol above and no ice and no radiation
git checkout run_scripts/westbrook_illingworth_nml.in
git checkout run_scripts/westbrook_illingworth_pamm_nml.in

sed -e "s/n_read(1,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(1,1:4)     = 250e6, 250e6, 250e6,0.1e6,/" run_scripts/westbrook_illingworth_pamm_nml.in > /tmp/namelist.tmp
sed -e "s/n_read(2,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(2,1:4)     = 250e6, 250e6, 250e6,0.1e6,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
sed -e "s/z_read(1:4)       = 0.,740,3260,4000/z_read(1:4)       = 0.,740,3700,3701/" /tmp/namelist.tmp2 > run_scripts/westbrook_illingworth_pamm_nml.in


sed -e "s/nm1%radiation=.true./nm1%radiation=.false./" run_scripts/westbrook_illingworth_pamm_nml.in > /tmp/namelist.tmp
sed -e "s/nm1%ice_flag=.true./nm1%ice_flag=.false./" /tmp/namelist.tmp /tmp/namelist.tmp2

cp /tmp/namelist.tmp2 run_scripts/westbrook_illingworth_pamm_nml.in


mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml.in

mv /tmp/output_pc01.nc /tmp/output0012.nc


git checkout run_scripts/westbrook_illingworth_nml.in
git checkout run_scripts/westbrook_illingworth_pamm_nml.in

