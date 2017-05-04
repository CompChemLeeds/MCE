#! /bin/bash

# Data combination script  
# Written by C. Symonds, 27/04/15

# This script combines the normpop files from multiple partial runs, requiring a list
# of the required folders WHICH SHOULD BE IN ORDER! The script calls an averaging 
# program <<subavrg.exe>> for static stepsizes, and the interpolation program 
# <<interpolate.exe>> to result in a final normpop.out file for each folder as would
# be needed by the collate.sh script (which is called at the end). The folderlist.dat
# file is built gradually by the restart.sh script

# Check that the folder list exists
if [[ ! -f folderlist.dat ]]; then
  echo "Folder list is missing! Aborting"
  exit 1
fi
runf=$PWD

IFS=$'\n' read -d '' -r -a lines < folderlist.dat # Read folder list into array
nelements=${#lines[@]}                            # Number of elements in array
max_index=$[$nelements-1]                         # Maximum array index

# Find the name of the temporary folder for combining data and creating it if needed
outf1=${lines[0]}
outf2=${outf1##*/}
cd ../
if [ -d $outf2 ]; then
  echo "The temp folder for concatenating the normpop files $outf2 already exists!"
  echo "Aborting."
  exit 1
else
  mkdir $outf2
  cd $outf2
  outf=$PWD
  cd ../run
fi

# Find the number of sub-folders in each folder
j=0
for i in ${lines[@]}; do
  cd $i
  folnum[j]=$( ls -d */ | wc -l )
  (( j++ ))
done

# Ensure that the numbers of sub-folders in each folder are the same
j=0
for i in ${folnum[@]}; do
  if [ $i -ne ${folnum[0]} ]; then
    echo "Error! Expected to find ${folnum[0]} subdirectories in ${lines[j]}."
    echo "Only found $i"
    exit 1
  fi
  (( j++ ))
done

folders=${folnum[0]}

cd ${lines[0]}/1-run
num1=$( grep "^Repeats" input.dat )
num2=${num1#* }
steptmp=$( grep "^step" prop.dat )
step=${steptmp#* }

# Find the end time of the simulation
endtmp=$( grep "^time_end" prop.dat )
end=${endtmp#* }
end2=$( echo $end | sed -e 's/d/e/g' ) # Change from, eg 1.0d0 to 1.0e0
end=$( echo $end2 | awk '{ printf "%#9.8E", $1 }' ) # Put in a common form

# Collect all the partial normpop files and put them in the temp folder
for i in `seq $folders`; do
  for k in `seq -f "%03g" $num2`; do
    p=1
    for j in ${lines[@]}; do
      cd $j/$i-run
      if [ -f normpop-$k.out ]; then
        cp normpop-$k.out $outf/normpop-${k}_${i}_${p}.out
        (( p++ ))
      fi
    done
  done
done

# Concatenate normpop files, omitting the headers
cd $outf
for k in `seq -f "%03g" $num2`; do
  for i in `seq $folders`; do
    for p in normpop-${k}_${i}_*.out; do
      q1=${p%.out}
      q=${q1##*_}
      if [ $q -eq 1 ]; then
        cp $p normpop-${k}_${i}.out
      else
        tail -n +5 -q $p >> normpop-${k}_${i}.out
      fi
    done
  done
done

# Checks to make sure normpop file ends at simulation end time
for k in `seq -f "%03g" $num2`; do
  for i in `seq $folders`; do
    tmp=$( tail -n 1 normpop-${k}_${i}.out )
    tmp2=${tmp#*  }
    tmp=${tmp2%% *}
    endtmp=$( echo $tmp | awk '{ printf "%#9.8E", $1 }' )
    if [[ $endtmp != $end ]]; then
      echo "The file normpop-${k}_${i}.out does not finish at the expected end time."
      echo "This file will be deleted"
      rm normpop-${k}_${i}.out
    fi
  done
done

# Move the combined normpop file to the final partial folder, and find number of
# columns in the normpop file (used as an argument for averaging)
n=0
path=${lines[max_index]}
cp -rf $path ${path}_bak
for k in `seq -f "%03g" $num2`; do
  for i in `seq $folders`; do
    if [ -f normpop-${k}_${i}.out ]; then
      n=$[$n+1]
      if [ $n -eq 1 ]; then
        cols=$( awk '{print NF}' normpop-${k}_${i}.out | sort -nu | tail -n 1 )
      fi
      mv normpop-${k}_${i}.out ${lines[max_index]}/$i-run/normpop-${k}.out
    fi
  done
done

cd $runf
rm -rf $outf    # Remove the temporary folder

# Assign the arguments for the averaging/interpolation programs
path=${lines[max_index]}
NUMBER=${path##*-}
reps=$[$num2*$folders]

# Run the interpolation/averaging program
for i in `seq $folders`; do
  if [ $step == "adaptive" ]; then
    cp interpolate.exe $path/$i-run
    cd $path/$i-run
    ./interpolate.exe $num2 $cols # Interpolates and averages data for subfolder
    cd $runf
  else
    cp subavrg.exe $path/$i-run
    cd $path/$i-run
    ./subavrg.exe $num2 $cols     # Averages normpop-*.out files in the subfolder
    line=$( cat normpop.out | wc -l )
    if [ $line -eq 3 ]; then rm normpop.out; fi
    cd $runf
  fi
done

# Create the results script and run it to combine the files from all sub-folders
echo "./collate.sh $path $reps $folders $NUMBER "'$0' > result.sh
chmod u+x result.sh
./result.sh
