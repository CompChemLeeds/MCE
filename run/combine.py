# Data combination script 
# Written by O. A. Bramley 04/01/21
#
# This is based on a bash cript written by C. Symonds. It combines the normpop files from
# multiple partial runs. The script calls an averaging program subavrg.exe for static 
# stepsizes and the interpolation program, interpolate.ex to produce a final normpop.out
# file for each folder as would be needed by the collate.py script which is called at the 
# end.

import shutil
import subprocess
import run
import glob
import os
import sys
import inputs

nodes=run.nodes
EXDIR1=os.getcwd()

for i in range(nodes):
    filedir=EXDIR1+"/run-"+str(i+1)
    counter=len(glob.glob1(filedir,"normpop-*.out"))
    if(i==0):
        counterprev=counter
    else:
        if(counter==counterprev):
            counterprev=counter
        else:
            sys.exit("Folder "+str(i+1)+" has "+str(counter)+" normpop files previous folders had "+str(counterprev))


time_end=inputs.prop['time_end']
size=len(time_end)
time_end=float(time_end[:size-3])*10**1
dtinit=inputs.prop['dtinit']
size=len(dtinit)
dtinit=float(dtinit[:-3])*10**-3
clone_freq=float(inputs.clone['clon_freq'])

if(inputs.clone["Cloning"]=='V1'):
    roar=2**(((time_end/dtinit)/clone_freq)-1)
    repeats=run.repeats*roar
else:
    repeats=run.repeats
    if(repeats!=(counter*nodes)):
        sys.exit("Expected "+str(repeats)+" repeats but got "+str(counter)+" per folder")

if os.path.exists(EXDIR1+"/combination_backup"):
    sys.exit("temporary combining folder already exists")
else:
     os.mkdir(EXDIR1+"/combination_backup")

tempDir=EXDIR1+"/combination_backup/"
for i in range(nodes):
    SUBDIR=EXDIR1+"/run-"+str(i+1)
    for j in range(counter):
        val=f"{(j+1):04d}"
        shutil.copy2(SUBDIR+"/normpop-"+val+".out",tempDir+"normpop-"+str(j+1)+"_"+str(i+1)+".out")

if(inputs.prop['step']=='adaptive'):
    prog='/interpolate.exe'
else:
    prog='/subavrg.exe'

for i in range(nodes):
    SUBDIR=EXDIR1+"/run-"+str(i+1)
    shutil.copy2(EXDIR1+prog,SUBDIR)

arg=str(int(repeats/nodes))
#print(counter)
for i in range(nodes):
    SUBDIR=EXDIR1+"/run-"+str(i+1)
    os.chdir(SUBDIR)
    subprocess.run(["."+prog,arg,"13"])
    os.chdir(EXDIR1)

#Builds result file
result=open(EXDIR1+"/result.sh","r+")
result_content=result.read()
dollar=0
apostrophe=0
for i in range(len(result_content)):
    if(result_content[i]=='$'):
        dollar=i
    if(result_content[i]=="'"):
        apostrophe=i
        break

new_string=result_content[0:dollar+5]+str(int(repeats))+" "+str(nodes)+result_content[apostrophe-1:]
result.seek(0)
result.truncate()
result.write(new_string)
result.close()
subprocess.run(['chmod', 'u+x', EXDIR1+'/result.sh'])
