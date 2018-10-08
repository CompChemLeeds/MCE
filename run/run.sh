#! /bin/bash
#####################################################################################
#
#    Bash Script for Parallel Open MP execution of the MCE / CCS program           
#    Written by C. Symonds                                                   21/03/17
#
#    This script has been built to be run using SGE. Use of another job submission
#    system may require modifications to certain parts of the script. The script is
#    designed to compile, copy all the required files into an execution folder, and
#    submit the program as a job. Various checks, output handling, parameter setting
#    and module loading procedures are also included.
#
#    To run, you must give three arguments at execution:
#       1) The number of repeats
#       2) The number of parallel threads per folder/node (max 8)
#       3) The number of folders/nodes
#
#    As there are module load commands, the script must be run as source if running
#    on ARC1/ARC2/POLARIS. This means that for a run of 128 repeats using 8 cores 
#    per node and 4 nodes, the execution command would be:
#
#    source ./run.sh 128 8 4 for ARC, or
#
#    ./run.sh 128 8 4 elsewhere
#
#    For portability the program is compiled using gfortran. The ifort 12.1 compiler 
#    results in problems due to a known bug with allocation through subroutines when
#    openMP is enabled, however the intel compiler can be used if not version 12.1
#    and this would require modifications to the makefile in the build folder.  
#    As such, any changes made to the program should take this into account and 
#    endevour to utilise intrinsic functions rather than intel specific or GNU 
#    specific extensions where possible. PGf95 is not supported for compilation.
#
#    Due to the nature of Shared Memory environments, the execution of the program 
#    is split into separate folders, each one containing a separate executable file 
#    which runs on a number of parallel threads defined in argument 2. Since on arc1 
#    there are 8 cores per node, this number of parallel threads is limited to a 
#    maximum of 8. Parallel execution is carried out through the repeat running of 
#    the basis set propagation. As such, argument 1 (the total number of repeats) 
#    must be equal to an integer multiple of 2 * threads * folders. The factor of 2 
#    originates from the propagation of a basis set followed by proagation of the 
#    complex conjugate of that basis set. If however conjugate propagation is not
#    needed and the conjugate repeats flag is disabled in the input.dat file, the 
#    limit will be set at threads * folders. If gnuplot is present on the computer 
#    then a set of graphs are created through the automatic creation and execution 
#    of gnuplot scripts. These scripts, and the data being plotted, are generated 
#    by the avrgpops.exe and timehist.exe executable files, and different plots can
#    be generated through changing the source file of the same name. Each running 
#    of this script runs the makefile to recompile the program if necessary. 
#
#    If you are running with a precalculated set of basis sets, put the output files 
#    from the basis set calculation in the run folder before execution and conjugate 
#    propagation should be disabled.
#
#####################################################################################

#This section checks the input arguments for errors.
grep -i "Conjugate_Repeats YES" input.dat > /dev/null
CRchk=$?
# Allow grid altering calibration test on chmlin45 (change hostname to disable)
if [[ -n $( echo $HOSTNAME | fgrep -e "chmlin451" ) ]]; then
  grdalt=1
else
  grdalt=0
fi
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
  echo "Not enough runs selected. Must be 1 or greater"
  exit 1
elif [[ $2 -lt 1 ]]; then
  echo "Not enough threads selected. Must be 1 or greater"
  exit 1
elif [[ $2 -gt 16 ]]; then
  echo "Too many threads selected. Maximum of 8 available"
  exit 1
elif [[ $3 -lt 1 ]]; then
  echo "Not enough folders selected. Must be 1 or greater"
  exit 1
elif [[ $3 -gt 100 ]]; then
  echo "Too many folders! Maximum of 100 simultaneous job submissions allowed!"
  exit 1
elif [[ $(( $1/$3 )) -ge 5000 ]]; then
  echo "Too many repeats per folder! Must be less than 5000!"
  exit 1
elif [[ $(( $1%(2*$2*$3) )) -ne 0 && $CRchk == 0 ]]; then
  echo "Number of repeats not valid for conjugate repetition"
  echo "Should be an integer multiple of 2 * cores * folders. Check input.dat"
  echo "CRchk = $CRchk and the factor is $[2*$2*$3]"
  exit 1
elif [[ $(( $1%($2*$3) )) -ne 0 ]]; then
  echo "Number of repeats must be an integer multiple of cores * folders"
  exit 1
elif [[ $(( $2*$3 )) -gt 100 ]]; then
  echo "Total number of cores should stay below 100"
  exit 1
else
  echo "Arguments checked"
  if [[ -n $( echo $HOSTNAME | fgrep -e "arc3" -e "polaris" -e "arc2" ) ]]; then
    HPCFLG=1
  else
    HPCFLG=0
  fi
  if [[ $HPCFLG -eq 1 ]]; then module load mkl; fi    # module loads for HPC systems
  echo "Running Starting"
fi

REPS=$(( $1/$3 ))         # number of repeats per folder
CORES=$2                  # number of threads required per folder
FOLDERS=$3                # number of folders
RUNF=$PWD                 # run folder
NUMBER=$RANDOM            # random number for unique identifier of current run
FILE="MCE$NUMBER.sh"      # name of the job submission script

# Checks for the build folder (will contain the makefiles)
if [[ -d ../build ]]; then
  cd ../build
  BUILD=$PWD
  cd $RUNF
else
  echo "No build folder! Aborting"
  exit 1
fi

#Check for existance of SGE on the machine
if [[ ! -z $( command -v qstat ) ]]; then HSTFLG=1; else HSTFLG=0; fi

#Determines the Execution directory
if [[ $HPCFLG -eq 0 ]]; then
  cd ..
  if [[ ! -d "EXEC" ]]; then
    mkdir EXEC
    cd EXEC
  else
    cd EXEC
  fi 
  EXDIR1=$PWD    # execution folder
else
  EXDIR1="/nobackup/$LOGNAME"
fi

folseq=( `seq 1 $FOLDERS` )    # The sequence of sub-folders, saved in an array

# Determine the correct makefile for the system and run make
cd $BUILD
if [[ HPCFLG -eq 1 ]]; then
  cp makefile_arc Makefile
  make
else
  cp makefile_chmlin Makefile
  make
fi 
cp *.exe $RUNF 
cd $RUNF
if [[ $? -ne 0 ]]; then 
  echo "Compilation Error! Exitting"
  exit 1
fi

cp input.dat input2.dat     # Created so it can be edited without conflicts

method=`grep -i "^method" input.dat`
if [[ $? != 0 ]]; then
  echo "Could not read the method from input.dat. Exitting"
  exit 1
fi
method=${method#* }
if [[ $method == "MCE12" ]]; then  # If method is MCE12, two runs are done
  k=2
else
  k=1
fi
methseq=( `seq 1 $k` )   # Array for number of executions (>1 element only for MCE12)

freqflg=`grep -i "^freqflg" input.dat`
if [[ $? != 0 ]]; then
  echo "Could not read the frequency flag from input.dat. Exitting"
  exit 1
fi
freqflg=${freqflg#* }
if [[ $freqflg != 0 && $freqflg != 1 ]]; then
  echo "Frequency flag is neither 1 or 0. Read a value of $freqflg"
  exit 1
fi

if [[ $freqflg -eq 1 && ! -f freq.dat ]]; then ./integrator.exe; fi

sed -i "s/^Repeats.*/Repeats $REPS/g" input2.dat   # Writes number of repeats per
grep "^Repeats $REPS" input2.dat > /dev/null       # folder to input.dat file
if [[ $? != 0 ]]; then
  echo "Could not change the number of repeats in input.dat. Exitting"
  exit 1
fi

# The name of the output folder is taken from the input.dat file
# If it is "default" then the folder name is made from the system, 
# the method and a random number, else it is made from the system,
# the method and whatever string is read from the line in input.dat.
# Determination of the folder name is carried out later, written to
# the variable $outfol2

outfol=`grep -i "^Runfolder" input.dat`
if [[ $? != 0 ]]; then
  echo "Could not read the execution folder from input.dat. Exitting"
  exit 1
fi
outfol=${outfol#* }
echo $outfol | grep -i "default" > /dev/null
outdef=$?

sys=`grep -i "^System:" input.dat`
if [[ $? != 0 ]]; then
  echo "Could not read the system from input.dat. Exitting"
  exit 1
fi
sys=${sys#* }

for a in "${methseq[@]}"; do

  # Determination of the output folder name
  if [[ $k == 2 ]]; then
    sed -i "s/^method.*/method MCEv$a/g" input2.dat
    if [[ $outdef == 0 ]]; then
      outfol2="MCEv$a-$sys-$NUMBER"
    else
      outfol2="MCEv$a-$sys-$outfol"
    fi
  else
    if [[ $outdef == 0 ]]; then
      outfol2="$method-$sys-$NUMBER"
    else
      outfol2="$method-$sys-$outfol"
    fi
  fi
  
  # Create output directory and all sub-folders
  EXDIR="$EXDIR1/$outfol2"
  if [[ ! -d $EXDIR ]]; then mkdir $EXDIR; fi 
  for i in "${folseq[@]}"; do
    SUBDIR="$EXDIR/$i-run"
    if [[ ! -d "$SUBDIR" ]]; then
      mkdir "$SUBDIR"
    else
      cd "$SUBDIR"
      if [[ "$(ls -A )" ]]; then rm *.*; fi      #remove old run files
    fi
  done
  
  # Check to see if basis set generation is enabled
  cd "$RUNF"
  chk=`grep "gen YES" input.dat`
  if [[ ! -z ${chk} ]]; then 
    gen=1
  else
    gen=0
  fi
  
  # Create the job submission file
  echo "#$ -cwd -V" > $FILE        # Run in CWD, to email at end of run add "-m e -M <email@address>"
  if [[ $CORES -ne 1 && $HPCFLG -eq 1 ]]; then 
    echo "#$ -pe smp $CORES" >> $FILE   # Use shared memory parallel environment
  fi
  if [[ $HPCFLG -eq 1 ]]; then
    echo "#$ -l h_rt=40:00:00" >> $FILE # Maximum allowed runtime (keep < 48h)
  fi
  echo "#$ -l h_vmem=4G" >> $FILE       # Allocated virtual memory
  echo "#$ -t 1-$FOLDERS" >> $FILE      # Set up job array
#  echo "#$ -tc 40" >> $FILE            # Maximum number of simultaneous running jobs
  echo "date" >> $FILE
  echo "cd $EXDIR/"'$SGE_TASK_ID'"-run/" >> $FILE
  echo "echo "'"Running on $HOSTNAME in folder $PWD"' >> $FILE 
  if [[ $HPCFLG -eq 1 ]]; then
    echo "module load mkl" >> $FILE     # Load mkl linear algebra set
  fi
  echo "time ./MCE.exe" >> $FILE             # Run program
  echo "date" >> $FILE
  
  if [[ $grdalt -eq 1 ]]; then ./gridchanger.sh; fi # Makes calibration input files
  
  for i in "${folseq[@]}"; do
    
    # Copy the input files / executables into each subfolder
    SUBDIR="$EXDIR/$i-run"
    cd "$RUNF"
    if [[ $grdalt -eq 1 ]]; then 
      cp ./calibinputs/input.$i $SUBDIR/input.dat
      cp inham.dat MCE.exe prop.dat $SUBDIR/
      if [[ $freqflg == 1 ]]; then 
        for x in `seq -f "%03g" 1 $REPS`; do
          cp freq.dat $SUBDIR/freq${x}.dat
        done
      fi
    else
      cp inham.dat input2.dat MCE.exe prop.dat $SUBDIR/
      if [[ $freqflg == 1 ]]; then 
        for x in `seq -f "%03g" 1 $REPS`; do
          cp freq.dat $SUBDIR/freq${x}.dat
        done
      fi
      mv $SUBDIR/input2.dat $SUBDIR/input.dat
    fi
    
    # Copy previously calculated basis set files if needed
    if [[ $gen -eq 0 ]]; then
      if [[ $method == "AIMC-MCE2" ]]; then
        if [[ -f "Outbs-001-00000-0_$i.out" || -f "Outbs-0001-00000-0_$i.out" ]]; then
          echo "Outbs-0001-00000-0_$i.out found in $PWD"
          for x in Outbs-*_$i.out; do
            cp $x $SUBDIR/${x%_$i.out}.out
          done
        else
          echo "Outbs-001-00000-0_$i.out not found in $PWD"
          echo "For AIMC-MCE second pass, all relevant input bases must be present"
          exit 1
        fi
        if [[ -f "Clonetrack-001_$i.out" || -f "Clonetrack-0001_$i.out" ]]; then
          echo "Clonetrack-001_$i.out found in $PWD"
          for x in Clonetrack-*_$i.out; do
            cp $x $SUBDIR/${x%_$i.out}.out
          done
        else
          echo "Clonetrack-001_$i.out not found in $PWD"
          echo "For AIMC-MCE second pass, the cloning tracking file must be present"
          exit 1
        fi 
      else
        if [[ -f "Outbs-001_$i.out" || -f "Outbs-0001_$i.out" ]]; then 
          echo "Outbs-0001_$i.out found in $PWD"
          for x in Outbs-*_$i.out; do
            cp $x $SUBDIR/${x%_$i.out}.out
          done
        else
          echo "Outbs-0001_$i.out not found in $PWD"
          echo "For propagation to occur without basis set generation,"
          echo "all relevant input bases must be present"
          exit 1
        fi
      fi
      
      # Copy cloning array files (used for restarted previous run
      clone=`grep -i "^Cloning" input.dat`
      if [[ $? != 0 ]]; then
        echo "Could not read the cloning flag from input.dat. Exitting"
        exit 1
      fi
      clone=${clone#* }
      if [[ $clone == "yes" ]]; then
        if [[ -f "clonearr-001_$i.out" || -f "clonearr-0001_$i.out" ]]; then 
          echo "clonearr-001_$i.out found in $PWD"
          for x in clonearr-*_$i.out; do
            cp $x $SUBDIR/${x%_$i.out}.out
          done
        else
          echo "clonearr-001_$i.out not found in $PWD"
          echo "For cloned propagation to occur without basis set generation, "
          echo "prior cloning information is needed"
          exit 1
        fi
      fi 
    fi
    
    # If not running on a SGE machine, run the program directly
    if [ $HSTFLG -eq 0 ]; then
      cd $SUBDIR/
      echo "Program Executing in $EXDIR"
      if [[ $CORES -ne 1 ]]; then export OMP_NUM_THREADS=$CORES; fi
      ./MCE.exe #&> $FILE.o1 &
      cd $RUNF
    fi
  done

  # If running on an SGE machine, submit the job array  
  mv $FILE $EXDIR
  echo "$RUNF/collate.sh $EXDIR $1 $3 $NUMBER "'$0' > $EXDIR/result.sh
  chmod u+x $EXDIR/result.sh
  cd $EXDIR
  if [[ $CORES -ne 1 ]]; then export OMP_NUM_THREADS=$CORES; fi
  if [[ $HSTFLG -eq 1 ]]; then
    qsub $FILE
  fi
  cd "$RUNF"
  
done

#Build the results.sh file to call the collate script, and clean up temp input file
echo "./collate.sh $EXDIR $1 $3 $NUMBER "'$0' > result.sh
rm input2.dat
chmod u+x result.sh
