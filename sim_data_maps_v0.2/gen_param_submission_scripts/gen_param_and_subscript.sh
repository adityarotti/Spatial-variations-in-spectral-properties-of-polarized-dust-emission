#!/bin/bash

datatype=(yr1 yr2)

#echo ${#whichsim[@]} # Returns the length of the array.

tot_nrlz=100
num_nodes=20
delta_num=$((tot_nrlz/num_nodes))


inparampath="../param/"

exec_name="run-gen-sim" 

pathout="../jobs/nn$num_nodes" ; mkdir $pathout
parampath=$pathout/param ; mkdir $parampath
wparampath=$pathout/which_param ; mkdir $wparampath
submitpath=$pathout/sub_script ; mkdir $submitpath

# Generates a set of parameter files given the number of nodes and the number of realization to be analyzed.
for ((i=0 ; i<=${#datatype[@]}-1; i++))
do
    for ((n=1 ; n<=num_nodes; n++))
    do
        nstart=$((1+delta_num*($n-1)))
        param_fname=$parampath/param_${datatype[$i]}_file${n}.in
        wparam_fname=$wparampath/param_${datatype[$i]}_file${n}
        echo "'"${param_fname#"."}"'" > $wparam_fname
        cp $inparampath/param_${datatype[$i]}.in $param_fname
        echo $nstart $delta_num "			--> Nstart, Nrlz" >> $param_fname
    done
done


# Generates the submission scripts to submit jobs. 
# Also submits the jobs on uncommenting the appropriate lines.
for ((i=0 ; i<=${#datatype[@]}-1; i++))
do
    for ((n=1 ; n<=num_nodes; n++))
    do
        wparam_fname=${wparampath#"."}/param_${datatype[$i]}_file${n}
        submit_fname=${submitpath}/submit_${datatype[$i]}_file${n}.j
        cp hopper_master_run_sim $submit_fname
        echo "#PBS -N "${datatype[$i]}_file${n} >> $submit_fname
	    echo "#PBS -e ./sysout/"${datatype[$i]}_file${n}".""$""PBS_JOBID.err" >> $submit_fname
        echo "#PBS -o ./sysout/"${datatype[$i]}_file${n}".""$""PBS_JOBID.out" >> $submit_fname
#        echo "export OMP_NUM_THREADS=8" >> $submit_fname
        echo "./"$exec_name " < " $wparam_fname >>  $submit_fname
    done
done
