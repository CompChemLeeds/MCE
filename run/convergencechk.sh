#! /bin/bash

# Script for the systematic testing of convergence through the modification of dimensionality and Nbf
# with the DL potential
# Written by C. Symonds, 21/04/15

# This script creates a set of input files for the MCE/CCS program
# by modifying the initial input file and changing the number of trajectories
# and/or the number of degrees of freedom. It then runs the run.sh script the required
# number of times to ensure that all versions are calculated.

# Inputs are the same as those desired FOR A SINGLE FOLDER

ndim_min=150
ndim_max=1500
nbf_min=150
nbf_max=1500
k=0

if [ ! -d calibinputs ]; then mkdir calibinputs; else cd calibinputs; rm *; cd ..; fi
if [ ! -f input2.dat ]; then cp input.dat input2.dat; fi

j=10
sed -i "s/^in_nbf.*/in_nbf $j/g" input2.dat
for ((i=ndim_min; i<=ndim_max; i=i+50)); do
 sed -i "s/^ndim.*/ndim $i/g" input2.dat   # Change the dimensionality
 sed -i "s/^Runfolder.*/Runfolder scaletst_${i}dim_${j}bf/g" input2.dat
 k=$[$k+1]
 cp input2.dat ./calibinputs/input.$k
done 

i=10
sed -i "s/^ndim.*/ndim $i/g" input2.dat   # Change the dimensionality
for ((j=nbf_min; j<=nbf_max; j=j+50)); do
 sed -i "s/^in_nbf.*/in_nbf $j/g" input2.dat
 sed -i "s/^Runfolder.*/Runfolder scaletst_${i}dim_${j}bf/g" input2.dat
 k=$[$k+1]
 cp input2.dat ./calibinputs/input.$k
done

rm input2.dat

echo "$k files created. Starting execution...."

for ((i=1; i<=$k; i++)); do 
  cp ./calibinputs/input.$i ./input.dat
  ./run.sh $1 $2 $3
done
