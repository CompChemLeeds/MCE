#$ -cwd -V
#$ -l h_vmem=4G
#$ -t 1-1
date
cd /home/ds/phy4cs/Dropbox/PhysChem/Github2017/MCE/EXEC/MCEv1-SB-testofamplitudes/$SGE_TASK_ID-run/
echo "Running on $HOSTNAME in folder $PWD"
time ./MCE.exe
date
