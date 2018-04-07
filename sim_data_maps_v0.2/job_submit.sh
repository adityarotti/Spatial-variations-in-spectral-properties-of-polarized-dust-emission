#!/bin/bash

datatype=(yr1 yr2)

#echo ${#datatype[@]} # Returns the length of the array.

tot_nrlz=100
num_nodes=20

pathout="./jobs"
submitpath=$pathout/nn${num_nodes}/sub_script

# Generates the submission scripts to submit jobs. 
# Also submits the jobs on uncommenting the appropriate lines.
prefix="../"
for ((i=0 ; i<=${#datatype[@]}-1; i++))
do
    for ((n=1 ; n<=num_nodes; n++))
    do
        qsub.serial ${submitpath}/submit_${datatype[$i]}_file${n}.j
        echo qsub.serial ${submitpath}/submit_${datatype[$i]}_file${n}.j
    done
done
