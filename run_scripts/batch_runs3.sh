#! /bin/bash

# batch runs



# run 1: high aerosol
git checkout run_scripts/westbrook_illingworth_nml_adiabatic.in
git checkout run_scripts/westbrook_illingworth_pamm_nml.in 
# change the humidity above
sed -e "s/nm1%rh_above=0.2,/nm1%rh_above=0.2,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
# change the jump in theta above
sed -e "s/nm1%th_jump=8.,/nm1%th_jump=5.,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
# change the gradient in theta above
sed -e "s/nm1%th_grad=0.0022,/nm1%th_grad=0.005,/" /tmp/namelist.tmp2 > /tmp/namelist.tmp3
cp /tmp/namelist.tmp3 run_scripts/westbrook_illingworth_nml_adiabatic.in


sed -e "s/n_read(1,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(1,1:4)     = 5000e6, 5000e6, 5000e6,5000e6,/" run_scripts/westbrook_illingworth_pamm_nml.in > /tmp/namelist.tmp
sed -e "s/n_read(2,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(2,1:4)     = 5000e6, 5000e6, 5000e6,5000e6,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
sed -e "s/n_read(3,1:4)     = 1e6, 1e6, 1e6,1e6,/n_read(3,1:4)     = 100e6, 100e6, 100e6,100e6,/" /tmp/namelist.tmp2 > run_scripts/westbrook_illingworth_pamm_nml.in 

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml_adiabatic.in 

mv /tmp/output_pc04.nc /tmp/output0001.nc

# run 2: lower aerosol
git checkout run_scripts/westbrook_illingworth_pamm_nml.in 
sed -e "s/n_read(1,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(1,1:4)     = 500e6, 500e6, 500e6,500e6,/" run_scripts/westbrook_illingworth_pamm_nml.in > /tmp/namelist.tmp
sed -e "s/n_read(2,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(2,1:4)     = 500e6, 500e6, 500e6,500e6,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
sed -e "s/n_read(3,1:4)     = 1e6, 1e6, 1e6,1e6,/n_read(3,1:4)     = 10e6, 10e6, 10e6,10e6,/" /tmp/namelist.tmp2 > run_scripts/westbrook_illingworth_pamm_nml.in 

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml_adiabatic.in 

mv /tmp/output_pc04.nc /tmp/output0002.nc

# run 3: lower aerosol
git checkout run_scripts/westbrook_illingworth_pamm_nml.in 

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml_adiabatic.in 

mv /tmp/output_pc04.nc /tmp/output0003.nc

# run 4: lower aerosol
git checkout run_scripts/westbrook_illingworth_pamm_nml.in 
sed -e "s/n_read(1,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(1,1:4)     = 125e6, 125e6, 125e6,125e6,/" run_scripts/westbrook_illingworth_pamm_nml.in > /tmp/namelist.tmp
sed -e "s/n_read(2,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(2,1:4)     = 125e6, 125e6, 125e6,125e6,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
sed -e "s/n_read(3,1:4)     = 1e6, 1e6, 1e6,1e6,/n_read(3,1:4)     = 0.1e6, 0.1e6, 0.1e6,0.1e6,/" /tmp/namelist.tmp2 > run_scripts/westbrook_illingworth_pamm_nml.in 

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml_adiabatic.in 

mv /tmp/output_pc04.nc /tmp/output0004.nc



# run 5: heymsfield and westbrook
sed -e "s/nm1%heyms_west=.true.,/nm1%heyms_west=.false.,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
cp /tmp/namelist.tmp run_scripts/westbrook_illingworth_nml_adiabatic.in 

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml_adiabatic.in 

mv /tmp/output_pc04.nc /tmp/output0005.nc


# run 6: no radiation
git checkout run_scripts/westbrook_illingworth_nml_adiabatic.in 
# change the humidity above
sed -e "s/nm1%rh_above=0.2,/nm1%rh_above=0.2,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
# change the jump in theta above
sed -e "s/nm1%th_jump=8.,/nm1%th_jump=5.,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
# change the gradient in theta above
sed -e "s/nm1%th_grad=0.0022,/nm1%th_grad=0.005,/" /tmp/namelist.tmp2 > /tmp/namelist.tmp3
cp /tmp/namelist.tmp3 run_scripts/westbrook_illingworth_nml_adiabatic.in

sed -e "s/nm1%radiation=.true.,/nm1%radiation=.false.,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
cp /tmp/namelist.tmp run_scripts/westbrook_illingworth_nml_adiabatic.in 

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml_adiabatic.in 

mv /tmp/output_pc04.nc /tmp/output0006.nc

# run 7: no divergence
git checkout run_scripts/westbrook_illingworth_nml_adiabatic.in 
# change the humidity above
sed -e "s/nm1%rh_above=0.2,/nm1%rh_above=0.2,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
# change the jump in theta above
sed -e "s/nm1%th_jump=8.,/nm1%th_jump=5.,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
# change the gradient in theta above
sed -e "s/nm1%th_grad=0.0022,/nm1%th_grad=0.005,/" /tmp/namelist.tmp2 > /tmp/namelist.tmp3
cp /tmp/namelist.tmp3 run_scripts/westbrook_illingworth_nml_adiabatic.in

sed -e "s/nm1%divergence=.true.,/nm1%divergence=.false.,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
cp /tmp/namelist.tmp run_scripts/westbrook_illingworth_nml_adiabatic.in 

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml_adiabatic.in 

mv /tmp/output_pc04.nc /tmp/output0007.nc





# run 8: mode 2 on
git checkout run_scripts/westbrook_illingworth_nml_adiabatic.in 
# change the humidity above
sed -e "s/nm1%rh_above=0.2,/nm1%rh_above=0.2,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
# change the jump in theta above
sed -e "s/nm1%th_jump=8.,/nm1%th_jump=5.,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
# change the gradient in theta above
sed -e "s/nm1%th_grad=0.0022,/nm1%th_grad=0.005,/" /tmp/namelist.tmp2 > /tmp/namelist.tmp3
cp /tmp/namelist.tmp3 run_scripts/westbrook_illingworth_nml_adiabatic.in

sed -e "s/nm1%mode2_ice_flag=0,/nm1%mode2_ice_flag=1,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
cp /tmp/namelist.tmp run_scripts/westbrook_illingworth_nml_adiabatic.in 

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml_adiabatic.in 

mv /tmp/output_pc04.nc /tmp/output0008.nc


# run 9: collisional breakup on
git checkout run_scripts/westbrook_illingworth_nml_adiabatic.in 
# change the humidity above
sed -e "s/nm1%rh_above=0.2,/nm1%rh_above=0.2,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
# change the jump in theta above
sed -e "s/nm1%th_jump=8.,/nm1%th_jump=5.,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
# change the gradient in theta above
sed -e "s/nm1%th_grad=0.0022,/nm1%th_grad=0.005,/" /tmp/namelist.tmp2 > /tmp/namelist.tmp3
cp /tmp/namelist.tmp3 run_scripts/westbrook_illingworth_nml_adiabatic.in

sed -e "s/nm1%coll_breakup_flag1=0,/nm1%coll_breakup_flag1=1,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
cp /tmp/namelist.tmp run_scripts/westbrook_illingworth_nml_adiabatic.in 

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml_adiabatic.in 

mv /tmp/output_pc04.nc /tmp/output0009.nc


# run 10: mode 2 and collisional breakup on
git checkout run_scripts/westbrook_illingworth_nml_adiabatic.in 
# change the humidity above
sed -e "s/nm1%rh_above=0.2,/nm1%rh_above=0.2,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
# change the jump in theta above
sed -e "s/nm1%th_jump=8.,/nm1%th_jump=5.,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
# change the gradient in theta above
sed -e "s/nm1%th_grad=0.0022,/nm1%th_grad=0.005,/" /tmp/namelist.tmp2 > /tmp/namelist.tmp3
cp /tmp/namelist.tmp3 run_scripts/westbrook_illingworth_nml_adiabatic.in

sed -e "s/nm1%mode2_ice_flag=0,/nm1%mode2_ice_flag=1,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
sed -e "s/nm1%coll_breakup_flag1=0,/nm1%coll_breakup_flag1=1,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
cp /tmp/namelist.tmp2 run_scripts/westbrook_illingworth_nml_adiabatic.in 

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml_adiabatic.in 

mv /tmp/output_pc04.nc /tmp/output0010.nc

# run 11: warm rain off
git checkout run_scripts/westbrook_illingworth_nml_adiabatic.in 
# change the humidity above
sed -e "s/nm1%rh_above=0.2,/nm1%rh_above=0.2,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
# change the jump in theta above
sed -e "s/nm1%th_jump=8.,/nm1%th_jump=5.,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
# change the gradient in theta above
sed -e "s/nm1%th_grad=0.0022,/nm1%th_grad=0.005,/" /tmp/namelist.tmp2 > /tmp/namelist.tmp3
cp /tmp/namelist.tmp3 run_scripts/westbrook_illingworth_nml_adiabatic.in

sed -e "s/nm1%wr_flag=.true.,/nm1%wr_flag=.false.,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
cp /tmp/namelist.tmp run_scripts/westbrook_illingworth_nml_adiabatic.in 

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml_adiabatic.in 

mv /tmp/output_pc04.nc /tmp/output0011.nc



# run 12: stochastic nucleation
git checkout run_scripts/westbrook_illingworth_nml_adiabatic.in 
# change the humidity above
sed -e "s/nm1%rh_above=0.2,/nm1%rh_above=0.2,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
# change the jump in theta above
sed -e "s/nm1%th_jump=8.,/nm1%th_jump=5.,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
# change the gradient in theta above
sed -e "s/nm1%th_grad=0.0022,/nm1%th_grad=0.005,/" /tmp/namelist.tmp2 > /tmp/namelist.tmp3
cp /tmp/namelist.tmp3 run_scripts/westbrook_illingworth_nml_adiabatic.in

sed -e "s/nm1%ice_nuc_flag=1,/nm1%ice_nuc_flag=2,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
sed -e "s/nm1%j_stochastic=0.5e-9,/nm1%j_stochastic=1e-9,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
cp /tmp/namelist.tmp2 run_scripts/westbrook_illingworth_nml_adiabatic.in 

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml_adiabatic.in 

mv /tmp/output_pc04.nc /tmp/output0012.nc

# run 13: stochastic nucleation
git checkout run_scripts/westbrook_illingworth_nml_adiabatic.in 
# change the humidity above
sed -e "s/nm1%rh_above=0.2,/nm1%rh_above=0.2,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
# change the jump in theta above
sed -e "s/nm1%th_jump=8.,/nm1%th_jump=5.,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
# change the gradient in theta above
sed -e "s/nm1%th_grad=0.0022,/nm1%th_grad=0.005,/" /tmp/namelist.tmp2 > /tmp/namelist.tmp3
cp /tmp/namelist.tmp3 run_scripts/westbrook_illingworth_nml_adiabatic.in

sed -e "s/nm1%ice_nuc_flag=1,/nm1%ice_nuc_flag=2,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
sed -e "s/nm1%j_stochastic=0.5e-9,/nm1%j_stochastic=1e-8,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
cp /tmp/namelist.tmp2 run_scripts/westbrook_illingworth_nml_adiabatic.in 

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml_adiabatic.in 

mv /tmp/output_pc04.nc /tmp/output0013.nc

# run 14: stochastic nucleation
git checkout run_scripts/westbrook_illingworth_nml_adiabatic.in 
# change the humidity above
sed -e "s/nm1%rh_above=0.2,/nm1%rh_above=0.2,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
# change the jump in theta above
sed -e "s/nm1%th_jump=8.,/nm1%th_jump=5.,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
# change the gradient in theta above
sed -e "s/nm1%th_grad=0.0022,/nm1%th_grad=0.005,/" /tmp/namelist.tmp2 > /tmp/namelist.tmp3
cp /tmp/namelist.tmp3 run_scripts/westbrook_illingworth_nml_adiabatic.in

sed -e "s/nm1%ice_nuc_flag=1,/nm1%ice_nuc_flag=2,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
sed -e "s/nm1%j_stochastic=0.5e-9,/nm1%j_stochastic=1e-7,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
cp /tmp/namelist.tmp2 run_scripts/westbrook_illingworth_nml_adiabatic.in 

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml_adiabatic.in 

mv /tmp/output_pc04.nc /tmp/output0014.nc





# run 15: shear 5
git checkout run_scripts/westbrook_illingworth_nml_adiabatic.in 
# change the humidity above
sed -e "s/nm1%rh_above=0.2,/nm1%rh_above=0.2,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
# change the jump in theta above
sed -e "s/nm1%th_jump=8.,/nm1%th_jump=5.,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
# change the gradient in theta above
sed -e "s/nm1%th_grad=0.0022,/nm1%th_grad=0.005,/" /tmp/namelist.tmp2 > /tmp/namelist.tmp3
cp /tmp/namelist.tmp3 run_scripts/westbrook_illingworth_nml_adiabatic.in

sed -e "s/nm1%param_vmax=2.,/nm1%param_vmax=5.,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
cp /tmp/namelist.tmp run_scripts/westbrook_illingworth_nml_adiabatic.in 

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml_adiabatic.in 

mv /tmp/output_pc04.nc /tmp/output0015.nc

# run 16: shear 10
git checkout run_scripts/westbrook_illingworth_nml_adiabatic.in 
# change the humidity above
sed -e "s/nm1%rh_above=0.2,/nm1%rh_above=0.2,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
# change the jump in theta above
sed -e "s/nm1%th_jump=8.,/nm1%th_jump=5.,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
# change the gradient in theta above
sed -e "s/nm1%th_grad=0.0022,/nm1%th_grad=0.005,/" /tmp/namelist.tmp2 > /tmp/namelist.tmp3
cp /tmp/namelist.tmp3 run_scripts/westbrook_illingworth_nml_adiabatic.in

sed -e "s/nm1%param_vmax=2.,/nm1%param_vmax=10.,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
cp /tmp/namelist.tmp run_scripts/westbrook_illingworth_nml_adiabatic.in 

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml_adiabatic.in 

mv /tmp/output_pc04.nc /tmp/output0016.nc

# run 17: shear 15
git checkout run_scripts/westbrook_illingworth_nml_adiabatic.in 
# change the humidity above
sed -e "s/nm1%rh_above=0.2,/nm1%rh_above=0.2,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
# change the jump in theta above
sed -e "s/nm1%th_jump=8.,/nm1%th_jump=5.,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
# change the gradient in theta above
sed -e "s/nm1%th_grad=0.0022,/nm1%th_grad=0.005,/" /tmp/namelist.tmp2 > /tmp/namelist.tmp3
cp /tmp/namelist.tmp3 run_scripts/westbrook_illingworth_nml_adiabatic.in

sed -e "s/nm1%param_vmax=2.,/nm1%param_vmax=15.,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
cp /tmp/namelist.tmp run_scripts/westbrook_illingworth_nml_adiabatic.in 

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml_adiabatic.in 

mv /tmp/output_pc04.nc /tmp/output0017.nc



# run 18: lawson process
git checkout run_scripts/westbrook_illingworth_nml_adiabatic.in 
# change the humidity above
sed -e "s/nm1%rh_above=0.2,/nm1%rh_above=0.2,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
# change the jump in theta above
sed -e "s/nm1%th_jump=8.,/nm1%th_jump=5.,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
# change the gradient in theta above
sed -e "s/nm1%th_grad=0.0022,/nm1%th_grad=0.005,/" /tmp/namelist.tmp2 > /tmp/namelist.tmp3
cp /tmp/namelist.tmp3 run_scripts/westbrook_illingworth_nml_adiabatic.in

sed -e "s/nm1%lawson=.false.,/nm1%lawson=.true.,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
cp /tmp/namelist.tmp run_scripts/westbrook_illingworth_nml_adiabatic.in 

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml_adiabatic.in 

mv /tmp/output_pc04.nc /tmp/output0018.nc


# run 19: no ice, no radiation
git checkout run_scripts/westbrook_illingworth_nml_adiabatic.in 
# change the humidity above
sed -e "s/nm1%rh_above=0.2,/nm1%rh_above=0.2,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
# change the jump in theta above
sed -e "s/nm1%th_jump=8.,/nm1%th_jump=5.,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
# change the gradient in theta above
sed -e "s/nm1%th_grad=0.0022,/nm1%th_grad=0.005,/" /tmp/namelist.tmp2 > /tmp/namelist.tmp3
cp /tmp/namelist.tmp3 run_scripts/westbrook_illingworth_nml_adiabatic.in

sed -e "s/nm1%ice_flag=.true.,/nm1%ice_flag=.false.,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
sed -e "s/nm1%radiation=.true.,/nm1%radiation=.false.,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
cp /tmp/namelist.tmp2 run_scripts/westbrook_illingworth_nml_adiabatic.in 

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml_adiabatic.in 

mv /tmp/output_pc04.nc /tmp/output0019.nc



# run 20: no aerosol above
git checkout run_scripts/westbrook_illingworth_nml_adiabatic.in
git checkout run_scripts/westbrook_illingworth_pamm_nml.in 
# change the humidity above
sed -e "s/nm1%rh_above=0.2,/nm1%rh_above=0.2,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
# change the jump in theta above
sed -e "s/nm1%th_jump=8.,/nm1%th_jump=5.,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
# change the gradient in theta above
sed -e "s/nm1%th_grad=0.0022,/nm1%th_grad=0.005,/" /tmp/namelist.tmp2 > /tmp/namelist.tmp3
cp /tmp/namelist.tmp3 run_scripts/westbrook_illingworth_nml_adiabatic.in



sed -e "s/n_read(1,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(1,1:4)     = 125e6, 125e6, 125e6,1e6,/" run_scripts/westbrook_illingworth_pamm_nml.in > /tmp/namelist.tmp
sed -e "s/n_read(2,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(2,1:4)     = 125e6, 125e6, 125e6,1e6,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
sed -e "s/z_read(1:4)       = 0.,740,3260,4000/z_read(1:4)       = 0.,740,3700,3701/" /tmp/namelist.tmp2 > /tmp/namelist.tmp3
cp /tmp/namelist.tmp3 run_scripts/westbrook_illingworth_pamm_nml.in

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml_adiabatic.in 

mv /tmp/output_pc04.nc /tmp/output0020.nc


# run 21: no recycling of aerosol
git checkout run_scripts/westbrook_illingworth_nml_adiabatic.in 
# change the humidity above
sed -e "s/nm1%rh_above=0.2,/nm1%rh_above=0.2,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
# change the jump in theta above
sed -e "s/nm1%th_jump=8.,/nm1%th_jump=5.,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
# change the gradient in theta above
sed -e "s/nm1%th_grad=0.0022,/nm1%th_grad=0.005,/" /tmp/namelist.tmp2 > /tmp/namelist.tmp3
cp /tmp/namelist.tmp3 run_scripts/westbrook_illingworth_nml_adiabatic.in

sed -e "s/nm1%recycle=.true.,/nm1%recycle=.false.,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
cp /tmp/namelist.tmp run_scripts/westbrook_illingworth_nml_adiabatic.in 

mpiexec -n 24 ./main.exe run_scripts/westbrook_illingworth_nml_adiabatic.in 

mv /tmp/output_pc04.nc /tmp/output0021.nc






