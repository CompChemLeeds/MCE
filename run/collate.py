#########################################################################################
# Data collection and collation script for the MCE/CCS program
# Written by O. A. Bramlye 20/10/2020

# This script is based off an earlier bash script written by C. Symonds
# This script collects all the data from the different subfolders, calls averaging
# programs to combine this data, then puts all relevant files in a results folder, 
# before deleting the raw data. By commenting out the marked sections of the script
# this raw data can be preserved. This script should be called only by use of the 
# results.sh file, which is automatically written by other scripts (either at the time
# of running, or when partial runs are combined using the combine.sh script) and so
# there should be no need to manually run this script alone unless the results.py
# script has been deleted/overwritten accidentally.
#########################################################################################

import sys
import os
import shutil
import glob
import subprocess

#Set purger to Y to delte raw data.
purger='N'

mcerunf=sys.argv[0]
path=sys.argv[1]
reps=sys.argv[2]
folders=sys.argv[3]
runfolder=sys.argv[4]
HPCFLG=sys.argv[5]
prop=sys.argv[6]

mcerunf=mcerunf[:-10]
os.chdir(mcerunf)
if not os.path.exists(path):
    sys.exit("Raw data values do not exist")

if os.path.exists("../"+runfolder):
    value=input("File already exists do you want to delete it? y/n\n")
    if(value=="y"):
        shutil.rmtree("../"+runfolder)
        os.mkdir("../"+runfolder)
    else:
        sys.exit("Results folder already exists. Delte/move it")
else:
    os.mkdir("../"+runfolder)
final_dir="../"+runfolder


for i in range(int(folders)):
    SUBDIR=path+"/run-"+str(i+1)
    print(SUBDIR)
    if not os.path.exists(SUBDIR):
        sys.exit("Error. Expected folder "+ SUBDIR+" does not exist")
    
    if os.path.exists(SUBDIR+"/normpop.out"):
        shutil.copy2(SUBDIR+"/normpop.out",final_dir+"/normpop_"+str(i+1)+".out")
    
    #Basis set files only retained when they have not been propagated
    if(prop=='NO'or prop=='no'or prop=='n'):
        for file in glob.glob(SUBDIR+'/Outbs-*.out'):
            shutil.copy(file,final_dir+"/"+file+"_"+str(i+1)+".out")

    if os.path.exists(SUBDIR+"/timehist.out"):
        shutil.copy2(SUBDIR+"/timehist.out",final_dir+"/timehist_"+str(i+1)+".out")
        shutil.copy2(SUBDIR+"/timesteps.out",final_dir+"/timesteps_"+str(i+1)+".out")

#Make short program to write out and present the run data
shutil.copy2(path+"/inputs.py",final_dir+"/inputs.py")
shutil.copy2(path+"/inham.py",final_dir+"/inham.py")

os.chdir(final_dir)
if os.path.exists("timehist_1.out"):
    with open(final_dir+"/timesteps.out",'wb') as wfd:
        for i in glob.glob(final_dir+"/timesteps_*.out"):
            with open(i,'rb') as fd:
                shutil.copyfileobj(fd,wfd)
    for file in glob.glob(final_dir+"/timesteps_*.out"):
        shutil.rmtree(file)
    
    shutil.copy2("../build/timehist.exe",final_dir+"/timehist.exe")
    timehist=subprocess.check_call(["./timehist.exe",str(folders),str(reps)])
    if(timehist==0):
        print("timehist.exe run correctly")
        os.remove("timehist.exe")
    else:
        print("timehist.exe did not run correctly")

if os.path.exists("normpop_1.out"):
    #with open("normpop_1.out","rb") as f:
    #   col=sum(1 for column in f)
    col=13
    shutil.copy2("../build/avrgpops.exe",final_dir+"/avrgpops.exe")
    avrgpops=subprocess.check_call(["./avrgpops.exe",str(folders),str(reps),str(col)])
    if(avrgpops==0):
        print("avrgpops has been run")
        os.remove("avrgpops.exe")
    else:
        print("avrgpops did not run correctly with arguements"+folders,reps,col)

gnu=subprocess.check_call(['which','gnuplot'])
if(gnu==0):
    for file in glob.glob("*.gpl"):
        subprocess.call(['gnuplot',file])

if(HPCFLG==1):
    for file in glob.glob(runfolder+"/*.oe*"):
        shutil.copy2(file,final_dir+"/"+file+".out")


purger=input("Do you want to delete raw data? y/n\n")
if(purger=='y'):
    shutil.rmtree(path)