#! /bin/bash

# Script for the systematic testing of trainj lengths and spacings
# Written for the testing of parameters for the spin boson potential
# Written by C. Symonds, 21/04/15

# This script creates a set of input files for the MCE/CCS program
# by modifying the initial input file and changing the size of the grid
# to accomodate different grid spacings such that the same area in phase space
# is covered. The number of basis functions is also modified so that it is correct.

# This script is currently only called on chmlin45, and only if the flag in run.sh 
# is enabled. This behaviour can be easily modified.

i=0

trsp=$1

if [ ! -d calibinputs ]; then mkdir calibinputs; else cd calibinputs; rm *; cd ..; fi

sed -i "s/^basis.*/basis SWARM/g" input.dat
sed -i "s/^method.*/method MCE12/g" input.dat
sed -i "s/^Cloning.*/Cloning no/" input.dat
sed -i "s/^max_cloning.*/max_cloning 4/" input.dat


#for trsp in 10; do
  cp input.dat input3.dat
  sed -i "s/^trainsp.*/trainsp $trsp/g" input3.dat   # Change the dimensions of the grid
  for defstp in 100; do
    cp input3.dat input4.dat
#    nbf=$[$defstp*11]
    nbf=100
    sed -i "s/^in_nbf.*/in_nbf $nbf/g" input4.dat
    sed -i "s/^def_stp.*/def_stp $defstp/g" input4.dat
    sed -i "/^Runfolder/ s/$/_${trsp}_clone/" input4.dat
    i=$[$i+1]
    mv input4.dat ./calibinputs/input.$i
  done
  rm input3.dat
#done 

echo "$i files created"
