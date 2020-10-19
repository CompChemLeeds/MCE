#$ -cwd -V 
#$ -pe smp 2 
#$ -l h_rt=40:00:00 
#$ -l h_vmem=4G 
#$ -t 1-1 
date 
cd ../EXEC/MCEv2-SB-t3/run-'$SGE_TASK_ID'/ 
echo "Running on $HOSTNAME in folder $PWD" 
module load mkl 
time ./MCE.exe 
date 
