#! /bin/bash

# Script for the systematic testing of grid spacings and compression parameters
# Written for the testing of parameters for the 3D Coulomb potential
# Written by C. Symonds, 21/04/15

# This script creates a set of input files for the MCE/CCS program
# by modifying the initial input file and changing the size of the grid
# to accomodate different grid spacings such that the same area in phase space
# is covered. The number of basis functions is also modified so that it is correct.

# This script is currently only called on chmlin45, and only if the flag in run.sh 
# is enabled. This behaviour can be easily modified.

qsize=140
psize=25
i=0

if [ ! -d calibinputs ]; then mkdir calibinputs; else cd calibinputs; rm *; cd ..; fi
if [ ! -f input2.dat ]; then cp input.dat input2.dat; fi

for grsp in 1.75; do
 dim1=$(perl -w -e "use POSIX; print ceil($qsize/$grsp), qq{\n}")     #Make integer
 dim2=$(perl -w -e "use POSIX; print ceil($psize/$grsp), qq{\n}")
 if [ $(( $dim1 % 2 )) == 1 ]; then dim1=$[$dim1+1]; fi  #Create even grid dimensions
 if [ $(( $dim2 % 2 )) == 1 ]; then dim2=$[$dim2+1]; fi
 dim3=$(( $dim1 * $dim2 ))
 sed -i "s/^qsizez.*/qsizez $dim1/g" input2.dat   # Change the dimensions of the grid
 sed -i "s/^psizez.*/psizez $dim2/g" input2.dat
 sed -i "s/^gridsp.*/gridsp $grsp/g" input2.dat
 sed -i "s/^in_nbf.*/in_nbf $dim3/g" input2.dat 
 for j in 500; do
  sed -i "s/^ALCMP.*/ALCMP $j/g" input2.dat
  sed -i "s/^Runfolder.*/Runfolder ${dim1}x${dim2}-$grsp-$j/g" input2.dat
   i=$[$i+1]
  cp input2.dat ./calibinputs/input.$i
 done
done 

echo "$i files created"
