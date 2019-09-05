#! /bin/bash
ARRAY1=(0.7 0.8 0.9 0.95 0.995) # rh above
ARRAY2=(0. 1. 5.) # theta jump above 


ELEMENTS1=${#ARRAY1[@]} # elements in first array
ELEMENTS2=${#ARRAY2[@]} # elements in second array

for (( i=0;i<$ELEMENTS1;i++)); do
	for (( j=0;j<$ELEMENTS2;j++)); do
	
        # Runs with the theta_flag switched on:
        echo ${ARRAY1[${i}]} ${ARRAY2[${j}]} ' theta_flag on'
        sed -e "s/rh_above=0.9/rh_above=${ARRAY1[${i}]}/" namelist_kelvin_helmholtz.in > /tmp/namelist.tmp
        sed -e "s/th_jump=0./th_jump=${ARRAY2[${j}]}/" /tmp/namelist.tmp > /tmp/namelist.run
        
        mpiexec -n 42 ../main.exe /tmp/namelist.run > std.out

        mv /tmp/output.nc /tmp/output_${i}_${j}_theta_on.nc
        
        
        
        # Runs with the theta_flag switched off:
        echo ${ARRAY1[${i}]} ${ARRAY2[${j}]} ' theta_flag off'
        sed -e "s/theta_flag=.true./theta_flag=.false./" /tmp/namelist.run > /tmp/namelist.tmp2
        sed -e "s/t_ctop=273./t_ctop=265./" /tmp/namelist.tmp2 > /tmp/namelist.run2
        
        mpiexec -n 42 ../main.exe /tmp/namelist.run2 > std.out

        mv /tmp/output.nc /tmp/output_${i}_${j}_theta_off.nc
 			
	done
done

