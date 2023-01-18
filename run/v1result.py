#########################################################################################
# Run folder data collection script for the MCEv1 cloning program
# Written by R. Brook 16/01/2023
# This script is based off collate.py and the result.sh file used to present results for general 
# MCE/CCS runs. This program should go through each of the run files and average over the different 
# normpop.out files that have already been condensed by the MCEv1 cloning. This reason makes it not 
# compatable for a run NOT USING MCEv1 cloning. After averaging the normpop file, a result file 
# should be created with the averaged file as well as graphs, the input files and gpl. This should be 
# run whenever cloning has occured instead of ./result.sh
#########################################################################################

import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
from pylab import cm
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

mcerunf=mcerunf[:-11]
os.chdir(mcerunf)
if not os.path.exists(path):
    sys.exit("Raw data values do not exist")

if os.path.exists("../"+runfolder):
    value=input("File already exists do you want to delete it? y/n\n")
    if(value=="y"or value=="Y"):
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
    fp=open('normpop_1.out')
    for count, line in enumerate(fp):
        pass
    fp.close()
    normp=np.zeros((count+1,13))
    for i in range(int(folders)):
        filename = final_dir+'/normpop_'+str(i+1)+'.out'
        print(filename)
        fp1=open(filename)
        for k, line in enumerate(fp1):
                if k>= 3:
                    list = line.split()
                    list = [float(p) for p in list]
                    list1=np.array(list)
                    normp[k,:]= np.add(normp[k,:],list1[:])
        fp1.close()
    
    normp = normp/int(folders)
    header = "Time Norm Re(ACF(t)) Im(ACF(t)) |ACF(t)| Re(Extra) Im(Extra) |Extra| Sum(HEhr) Pop1 Pop2 Pop1+Pop2 Pop2-Pop1 \n \n"
    np.savetxt(final_dir+"/normpop.out",normp[3:,:], delimiter=' ', header=  header, comments= '')
    mpl.rcParams['font.family']='Avenir'
    plt.rcParams['font.size']=18
    plt.rcParams['axes.linewidth']=2
    colors =cm.get_cmap('Set1',3)
    fig=plt.figure(figsize=(3.37,5.055))
    ax=fig.add_axes([0,0,2,1])
    x=normp[3:,0]
    y=normp[3:,12]
    ax.plot(x,y,linewidth=2, color=colors(0))
    ax.set_xlim(0,10)
    ax.set_ylim(-1,1)
    ax.set_ylabel('Population difference',labelpad=10)
    ax.set_xlabel('Time',labelpad=10)
    plt.savefig('popsdifftot.png',dpi=300, transparent=False,bbox_inches='tight') 
      
    fig2=plt.figure(figsize=(3.37,5.055))
    ax2=fig2.add_axes([0,0,2,1])
    x=normp[3:,0]
    y=normp[3:,10]
    ax2.plot(x,y,linewidth=2, color=colors(0))
    x=normp[3:,0]
    y=normp[3:,9]
    ax2.plot(x,y,linewidth=2, color=colors(1))
    ax2.set_xlim(0,10)
    ax2.set_ylim(0,1)
    ax2.set_ylabel('Population difference',labelpad=10)
    ax2.set_xlabel('Time',labelpad=10)
    plt.savefig('popstot.png',dpi=300, transparent=False,bbox_inches='tight') 





if(HPCFLG==1):
    for file in glob.glob(runfolder+"/*.oe*"):
        shutil.copy2(file,final_dir+"/"+file+".out")


purger=input("Do you want to delete raw data? y/n\n")
if(purger=='y'):
    shutil.rmtree(path)