#! /bin/bash

git checkout run_scripts/westbrook_illingworth_nml_adiabatic.in
git checkout run_scripts/westbrook_illingworth_pamm_nml.in

#sed -e "s/nm1%ice_nuc_flag=1,/nm1%ice_nuc_flag=2,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
#sed -e "s/nm1%j_stochastic=0.5e-9,/nm1%j_stochastic=1e-7,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
#cp /tmp/namelist.tmp2 run_scripts/westbrook_illingworth_nml_adiabatic.in

#sed -e "s/nm1%ice_flag=.true.,/nm1%ice_flag=.false.,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
#sed -e "s/nm1%radiation=.true.,/nm1%radiation=.false.,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
#cp /tmp/namelist.tmp2 run_scripts/westbrook_illingworth_nml_adiabatic.in

#sed -e "s/nm1%runtime= 28800./nm1%runtime= 86400./" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
#cp /tmp/namelist.tmp run_scripts/westbrook_illingworth_nml_adiabatic.in

#sed -e "s/output_pc04.nc/output_pc05.nc/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
#cp /tmp/namelist.tmp run_scripts/westbrook_illingworth_nml_adiabatic.in

# change the humidity above
sed -e "s/nm1%rh_above=0.2,/nm1%rh_above=0.2,/" run_scripts/westbrook_illingworth_nml_adiabatic.in > /tmp/namelist.tmp
# change the jump in theta above
sed -e "s/nm1%th_jump=8.,/nm1%th_jump=5.,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
# change the gradient in theta above
sed -e "s/nm1%th_grad=0.0022,/nm1%th_grad=0.005,/" /tmp/namelist.tmp2 > /tmp/namelist.tmp3


cp /tmp/namelist.tmp3 run_scripts/westbrook_illingworth_nml_adiabatic.in




# lower aerosol
git checkout run_scripts/westbrook_illingworth_pamm_nml.in
sed -e "s/n_read(1,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(1,1:4)     = 125e6, 125e6, 125e6,125e6,/" run_scripts/westbrook_illingworth_pamm_nml.in > /tmp/namelist.tmp
sed -e "s/n_read(2,1:4)     = 250e6, 250e6, 250e6,250e6,/n_read(2,1:4)     = 125e6, 125e6, 125e6,125e6,/" /tmp/namelist.tmp > /tmp/namelist.tmp2
sed -e "s/z_read(1:4)       = 0.,740,3260,4000/z_read(1:4)       = 0.,740,3700,3701/" /tmp/namelist.tmp2 > /tmp/namelist.tmp3
#sed -e "s/d_read(2,1:4)     = 150.e-9,150.e-9,150.e-9,150.e-9,/d_read(2,1:4)     = 250.e-9,250.e-9,250.e-9,250.e-9,/" /tmp/namelist.tmp2 > /tmp/namelist.tmp3
sed -e "s/n_read(3,1:4)     = 1e6, 1e6, 1e6,1e6,/n_read(3,1:4)     = 0.1e6, 0.1e6, 0.1e6,0.1e6,/" /tmp/namelist.tmp3 > run_scripts/westbrook_illingworth_pamm_nml.in

mpiexec -n 12 ./main.exe run_scripts/westbrook_illingworth_nml_adiabatic.in



git checkout run_scripts/westbrook_illingworth_nml_adiabatic.in
git checkout run_scripts/westbrook_illingworth_pamm_nml.in


