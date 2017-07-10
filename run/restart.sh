#! /bin/bash

# Simulation Restarting Script
# Written by C. Symonds, 27/04/15

# This script restarts a previous propagation. To use, call from last propagation set
# folder (eg in /nobackup folder). You should provide the correct number of repeats, 
# cores and folders which should match the parameters from the previous run. If 
# multiple simulations are being run, ensure that the correct folderlist.dat file and
# Outbsbackup directory are using the default names to ensure that there is no 
# overwriting of basis files or confusing of folder lists.

# Checking the input arguments

if [[ $# -ne 3 ]]; then
  echo "Error: Incorrect number of arguments! Arguments should be:"
  echo "       1: The total number of repeats"
  echo "       2: The number of threads required and "
  echo "       3: The number of folders"
  exit 1
elif [[ "`echo $1 | egrep ^[[:digit:]]+$`" = "" ]]; then  # checks for integer
  echo "Number of runs must be an integer"
  exit 1
elif [[ "`echo $2 | egrep ^[[:digit:]]+$`" = "" ]]; then  # checks for integer
  echo "Number of nodes must be an integer"
  exit 1
elif [[ "`echo $3 | egrep ^[[:digit:]]+$`" = "" ]]; then  # checks for integer
  echo "Number of folders must be an integer"
  exit 1
elif [[ $1 -lt 1 ]]; then
  echo "Not enough runs selected. Must be 2 or greater"
  exit 1
elif [[ $2 -lt 1 ]]; then
  echo "Not enough threads selected. Must be 1 or greater"
  exit 1
elif [[ $2 -gt 8 ]]; then
  echo "Too many threads selected. Maximum of 8 available"
  exit 1
elif [[ $3 -lt 1 ]]; then
  echo "Not enough folders selected. Must be 1 or greater"
  exit 1
elif [[ $(( $1/$3 )) -ge 500 ]]; then
  echo "Too many repeats per folder! Must be less than 500!"
  exit 1
elif [[ $(( $1%($2*$3) )) -ne 0 ]]; then
  echo "Number of repeats must be an integer multiple of cores * folders"
  exit 1
else
  folder=$PWD
  exec=${folder%/*}
  runf=${0%/restart.sh}                 
  folnum=$( ls -lR | grep ^d | wc -l )  # Find the number of folders in the folder
  
  # Copying the restart files into the run folder
  for i in `seq $folnum`; do
    cd ${i}-run
    if [ $i -eq 1 ]; then cp input.dat inham.dat prop.dat $runf; fi
    for j in Outbs-*; do
      cp $j $runf/${j%.out}_${i}.out
    done
    for j in clonearr-*; do
      cp $j $runf/${j%.out}_${i}.out
    done
    cd ..
  done
  
  cd $runf
  echo $folder >> folderlist.dat      # the folderlist is used for combining data
  restrtnum=$( cat folderlist.dat | wc -l )   # Number of entries in the folder list
  
  # Copy the restart files into a backup folder
  if [ ! -d $exec/Outbsbackup ]; then
    mkdir $exec/Outbsbackup
  fi
  for i in Outbs-*.out; do 
    cp $i /nobackup/phy4cs/Outbsbackup/${i%.out}_${restrtnum}.out
  done
  for i in clonearr-*.out; do 
    cp $i /nobackup/phy4cs/Outbsbackup/${i%.out}_${restrtnum}.out
  done
  
  # Ensure that the basis set is not recalculated
  chk=`grep "gen YES" input.dat`
  if [[ ! -z ${chk} ]]; then 
    sed -i "s/^gen.*/gen NO/g" input.dat
  fi
  
  source run.sh $1 $2 $3   # Restart the simulation
  rm Outbs* clonearr*      # Remove the restart files
fi
