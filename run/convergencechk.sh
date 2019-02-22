#! /bin/bash

# Script for the systematic testing of convergence through the modification of dimensionality and Nbf
# with the DL potential
# Written by C. Symonds, 21/04/15

# This script creates a set of input files for the MCE/CCS program
# by modifying the initial input file and changing the number of trajectories
# and/or the number of degrees of freedom. It then runs the run.sh script the required
# number of times to ensure that all versions are calculated.

# Inputs are the same as those desired FOR A SINGLE FOLDER

ndim_min=20
ndim_max=40
nbf_min=100
nbf_max=500
k=0

if [ ! -d calibinputs ]; then mkdir calibinputs; else cd calibinputs; rm *; cd ..; fi
if [ ! -f input2.dat ]; then cp input.dat input2.dat; fi

for i in 0 1; do
 sed -i "s/^qss.*/qss $i/g" input2.dat   # Change the qss flag
 for j in 100 200 400; do
  sed -i "s/^in_nbf.*/in_nbf $j/g" input2.dat
  sed -i "s/^Runfolder.*/Runfolder UBO_qss${i}_60dim_${j}bf/g" input2.dat
   k=$[$k+1]
  cp input2.dat ./calibinputs/input.$k
 done
done 

rm input2.dat

echo "$k files created. Starting execution...."

for ((i=1; i<=$k; i++)); do 
  cp ./calibinputs/input.$i ./input.dat
  ./run.sh $1 $2 $3
  mv ./result.sh result_${i}.sh
done
