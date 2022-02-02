#########################################################################################
#
#   Python Run script for Parallel Open MP execution of the MCE / CCS program
#   Written by O.A. Bramley                                     07/10/2020
#
#   This script is based of a similar script witted by C. Symonds using bash. This script
#   aims to simplify the running process and make the program more useable as python is
#   widely understood and should make modifications easier to implement.
#   The script is designed to compile, copy all reaquired files into an execution folder, 
#   and submit the program as a job. Included are various checks, output handling,  
#   parameter setting and module loading porcedures. This script  can also be used for 
#   for restarting a timed-out simulation by setting the restart paramter to 'YES'.
#   
#   To run the program variables must be set/checked in inputs.py, inham.py
#   The following arguemtns then have ot be set in the run folder
#       1) The number of repeats
#       2) The number of folders/nodes
#       3) The number of parallel cores per folder/node (max8)
#   
#   The propagation and basis set generation flags are also set in this file.
#   To restart a run open the run script in the execution folder
#   change the resart paramter to 'Yes' and run the copy of this script from the 
#   execution file. Starting a run with precalcualted basis funcions is not yet
#   possible for MCE12. 
#
#
#########################################################################################
import sys
import socket
import os
import subprocess
import getpass
import random
import shutil
import glob
import csv
import inham
import inputs

#########################################################################################
#                              VARIABLES TO SET FOR SIMULATION                          #
#########################################################################################

# Number of repeats 
repeats=1
# Number of nodes/folders
nodes=1
#Number of parallel cores per folder/node (max 8)
cores=1
# Name of running folder 
# Default : <method>-<system>-<random number> ie CCS-HP-31254
# Otherwise:  <method>-<system>-<runfolder string>
Runfolder='1rep'
# Generate Basis Set? YES/NO
gen='YES'
# Propagate Basis Set? YES/NO
prop='YES'
# Restart? YES/NO
# To restart a timedout run set to yes and rerun this script from the execution folder
restart='NO'
# Seed value for doing the random number routine- if you do not specify it 
# (leave default value of 0) will automatically generate one
SEED=0
#########################################################################################
#                                   END OF INPUTS                                       #
#########################################################################################
#                * NO NEED TO SCROLL FURTHER IF USING AS BLACKBOX *                     #
#########################################################################################


if __name__=="__main__":
    #Check basic arguements
    if(isinstance(repeats,int)==False):
        sys.exit("Number of repeats must be an integer")
    elif(isinstance(nodes,int)==False):
        sys.exit("Number of folders must be an integer")
    elif(isinstance(cores,int)==False):
        sys.exit("Number of parallel cores must be an integer")
    elif(repeats<1):
        sys.exit("Not enough runs selected. Must be 1 or greater")
    elif(nodes<1):
        sys.exit("Not enough nodes selected. Must be 1 or greater")
    elif(nodes>100):
        sys.exit("Too many nodes. Maximum of 100 simultaneous submisions")
    elif(cores>8):
        sys.exit("Too many cores selected. Maximum of 8 available")
    elif(cores<1):
        sys.exit("Not enough cores selected. Must be 1 or greater")
    elif((repeats/nodes)>5000):
        sys.exit("Too many repeats per folder. Must be less than 500")

    if(restart=="NO"):
        if((inputs.Conjugate_Repeats=='YES')and((repeats%(2*nodes*cores))!=0)):
            sys.exit("Number of repeats not valid for conjugate repetition. Should be integer multiple of 2*cores*nodes")
        elif((repeats%(nodes*cores))!=0):
            sys.exit("Number of repeats must be an integer multiple of cores*folders")
        elif(nodes*cores>100):
            sys.exit("Total number of cores should stay below 100")
        elif(inputs.systems['freqflg']not in{0,1}):
            sys.exit("Frequency flag msut be zero or 1")
        else:
            print("Arguments checked")
            Hostname=socket.gethostname()
            if(Hostname==("login1.arc4.leeds.ac.uk")):
                HPCFLG=1
            else:
                HPCFLG=0

        #Might need grid altering calibration test for chlin451 bash code
        #if [[ -n $( echo $HOSTNAME | fgrep -e "chmlin451" ) ]]; then
        #grdalt=1
        #else
        #grdalt=0
        #fi

        #Makes execution folder and run folder
        if(HPCFLG==0):
            if not os.path.exists("../EXEC"):
                os.mkdir("../EXEC")
            EXDIR="../EXEC"
        else:
            # subprocess.run(['module','load','mkl'])
            os.environ['LOGNAME']
            EXDIR="/nobackup/"+getpass.getuser()

        if(Runfolder=="Default"):
            Runfolder=inputs.method+"-"+inputs.systems["System"]+"-"+str(repeats)+"-"+str(nodes)+"-"+str(cores)
        else:
            Runfolder=inputs.method+"-"+inputs.systems["System"]+"-"+Runfolder


        if os.path.exists(EXDIR+"/"+Runfolder):
            value=input("File already exists do you want to delete it? y/n\n")
            if(value=='y'):
                shutil.rmtree(EXDIR+"/"+Runfolder)
            else:
                sys.exit("Runfolder already exists. Change the Runfolder name or delte/move it")
        
        os.mkdir(EXDIR+"/"+Runfolder)
        
        EXDIR1=EXDIR+"/"+Runfolder  

        mcerunf=os.getcwd()
        #Builds result file
        result=open(EXDIR1+"/result.sh","w")
        result.write("python "+mcerunf+"/collate.py $PWD "+(str(repeats))+" "+str(nodes)+" '"+Runfolder+"' "+(str(HPCFLG))+" '"+prop+"'")
        result.close()
        subprocess.run(['chmod', 'u+x', EXDIR1+'/result.sh'])

        #Copies input files
        shutil.copy2("inham.py",EXDIR1)
        shutil.copy2("inputs.py",EXDIR1)
        shutil.copy2("run.py",EXDIR1)
        shutil.copy2("combine.py",EXDIR1)

        
        #Makes the program input file
        if(inputs.method=="MCE12"):
            for i in range(2):
                with open('rundata'+str(i+1)+'.csv','w',newline='')as file:
                    writer = csv.writer(file)
                    writer.writerow([gen,prop,restart,inputs.cmprss,('MCEv'+str(i+1)),int(repeats/nodes),inputs.Conjugate_Repeats])
                    writer.writerow(inputs.systems.values())
                    writer.writerow(inputs.parameters.values())
                    writer.writerow(inputs.Train.values())
                    writer.writerow(inputs.clone.values())
                    writer.writerow(inputs.paramz.values())
                    writer.writerow(inham.EL.values())
                    writer.writerow(inputs.prop.values())
                    if(inputs.systems['System']=='MP'):
                        writer.writerow(inham.MP.values())
                    elif(inputs.systems['System']=='HP'):
                        writer.writerow(inham.HP.values())
                    else:
                        writer.writerow(inham.SB.values())
                shutil.copy2('rundata'+str(i+1)+'.csv',EXDIR1)  
        else:
            with open('rundata.csv','w',newline='')as file:
                writer = csv.writer(file)
                writer.writerow([gen,prop,restart,inputs.cmprss,inputs.method,int(repeats/nodes),inputs.Conjugate_Repeats])
                writer.writerow(inputs.systems.values())
                writer.writerow(inputs.parameters.values())
                writer.writerow(inputs.Train.values())
                writer.writerow(inputs.clone.values())
                writer.writerow(inputs.paramz.values())
                writer.writerow(inham.EL.values())
                writer.writerow(inputs.prop.values())
                if(inputs.systems['System']=='MP'):
                    writer.writerow(inham.MP.values())
                elif(inputs.systems['System']=='HP'):
                    writer.writerow(inham.HP.values())
                else:
                    writer.writerow(inham.SB.values())
            shutil.copy2("rundata.csv",EXDIR1)
        
        for file in glob.glob(mcerunf+"/*.csv"):
            os.remove(file)

        #Makes subfolders
        if(inputs.method=="MCE12"):
            os.mkdir(EXDIR1+"/MCEv1")
            os.mkdir(EXDIR1+"/MCEv2")
            for j in range(2):
                for i in range (nodes):
                    os.mkdir(EXDIR1+"/MCEv"+str(j+1)+"/run-"+str(i+1))         
        else:
            for i in range(nodes):
                path=os.path.join(EXDIR1,"run-"+str(i+1))
                os.mkdir(EXDIR1+"/run-"+str(i+1))
            

        #Selects the right make file and executes
        os.chdir("../build")
        if(HPCFLG==1):
            shutil.copy2("../build/makefile_arc","../build/Makefile")
            subprocess.run(["make"])
        else:
            shutil.copy2("../build/makefile_chmlin","../build/Makefile")
            subprocess.run(["make"])
        shutil.copy2("MCE.exe",EXDIR1)
        shutil.copy2("interpolate.exe",EXDIR1)
        shutil.copy2("subavrg.exe",EXDIR1)

        if(inputs.systems['freqflg']==1):
            if os.path.exists(mcerunf+"/freq.dat"):
                shutil.copy2(mcerunf+"/freq.dat",EXDIR1)
            else:
                subprocess.run(["./integrator.exe"])
                shutil.copy2("freq.dat",EXDIR1)

        os.chdir(EXDIR1)
        EXDIR1=os.getcwd()

        if(gen=='NO'): 
            if(inputs.method=="AIMC-MCE2"):
                if (glob.glob(mcerunf+"Outbs-001-00000-0_*.out")) or (glob.glob(mcerunf+"/Outbs-0001-00000-0_*.out")):
                    print("Outbs-0001-00000-0_*.out found in"+mcerunf)
                    for file in glob.glob(mcerunf+"/Outbs-*.out"):
                        shutil.copy2(file,EXDIR1)
                else:
                    sys.exit("Outbs-001-00000-0_*.out not found in runfolder For AIMC-MCE second pass, all relevant input bases must be present")
                if (glob.glob(mcerunf+"Clonetrack-001_*.out")) or (glob.glob(mcerunf+"/Clonetrack-0001_*.out")):
                    print("Clonetrack-0001_*.out found in"+mcerunf)
                    for file in glob.glob(mcerunf+"/Clonetrack-*.out"):
                        shutil.copy2(file,EXDIR1)
                else:
                    sys.exit("Clonetrack-001_*.out not found in runfolder For AIMC-MCE second pass, all relevant input bases must be present")
            else:
                if (glob.glob(mcerunf+"Outbs-001_*.out")) or (glob.glob(mcerunf+"/Outbs-0001_*.out")):
                    print("Outbs-0001_*.out found in"+mcerunf)
                    for file in glob.glob(mcerunf+"/Outbs-*.out"):
                        shutil.copy2(file,EXDIR1)
                else:
                    sys.exit("Outbs-001_*.out not found in runfolder.")
            if(inputs.clone['cloning']=='yes'):
                if (glob.glob(mcerunf+"Clonetrack-001_*.out")) or (glob.glob(mcerunf+"/Clonetrack-0001_*.out")):
                    print("Clonetrack-0001_*.out found in"+mcerunf)
                    for file in glob.glob(mcerunf+"/Clonetrack-*.out"):
                        shutil.copy2(file,EXDIR1)
                else:
                    sys.exit("Clonetrack-001_*.out not found in runfolder")

        if(inputs.method=="MCE12"):
            for j in range(2):
                for i in range (nodes):
                    shutil.copy2("MCE.exe","MCEv"+str(j+1)+"/run-"+str(i+1))
                    shutil.copy2('rundata'+str(j+1)+'.csv',"MCEv"+str(j+1)+"/run-"+str(i+1)+"/rundata.csv")
                    if(inputs.systems['freqflg']==1):
                        for k in range(repeats/nodes):
                            shutil.copy2("freq.dat","MCEv"+str(j+1)+"/run-"+str(i+1)+"/freq"+str(k+1)+".dat")
        else:
            for i in range (nodes):
                shutil.copy2("MCE.exe","run-"+str(i+1))
                shutil.copy2("rundata.csv","run-"+str(i+1))
                if(inputs.systems['freqflg']==1):
                    for k in range(repeats/nodes):
                        shutil.copy2("freq.dat","/run-"+str(i+1)+"/freq"+str(k+1)+".dat")
                if(gen=='NO'):
                    for file in glob.glob("Outbs-*.out"):
                        shutil.copy2(file,"/run-"+str(i+1))
                    if(inputs.clone['cloning']=='yes'):
                        for file in glob.glob("Clonetrack-*.out"):
                            shutil.copy2(file,"/run-"+str(i+1))
        


    elif(restart=='YES'):
        Hostname=socket.gethostname()
        if(Hostname==("login2.arc4.leeds.ac.uk")):
                HPCFLG=1
        else:
            HPCFLG=0

        EXDIR1=os.getcwd()

        if not os.path.exists("../Outbsbackup"):
                os.mkdir("../Outbsbackup")

        if(inputs.method=="MCE12"):
            for j in range(2):
                with open('rundata'+str(i+1)+'.csv','w',newline='')as file:
                    writer = csv.writer(file)
                    writer.writerow(['NO',prop,restart,inputs.cmprss,('MCEv'+str(j+1)),int(repeats/nodes),'NO'])
                    writer.writerow(inputs.systems.values())
                    writer.writerow(inputs.parameters.values())
                    writer.writerow(inputs.Train.values())
                    writer.writerow(inputs.clone.values())
                    writer.writerow(inputs.paramz.values())
                    writer.writerow(inham.EL.values())
                    writer.writerow(inputs.prop.values())
                    if(inputs.systems['System']=='MP'):
                        writer.writerow(inham.MP.values())
                    elif(inputs.systems['System']=='HP'):
                        writer.writerow(inham.HP.values())
                    else:
                        writer.writerow(inham.SB.values())
                for i in range(nodes):
                    p=1
                    q=1
                    for file in glob.glob('MCEv'+str(j+1)+'/run-'+str(i+1)+'/Outbs-*.out'):
                        shutil.copy2(file,'../Outbsbackup/Outbs-'+str(j+1)+'_'+str(p)+'.out_'+str(i+1)+'.out')
                        p=p+1
                    if((inputs.clone['clone'])!='no'):
                        for file in glob.glob('MCEv'+str(j+1)+'/run-'+str(i+1)+'/clonearr-*.out'):
                            shutil.copy2(file,'../Outbsbackup/clonearr-'+str(j+1)+'_'+str(p)+'.out_'+str(i+1)+'.out')
                            q=q+1   
                    shutil.copy2('rundata'+str(j+1)+'.csv',"MCEv"+str(j+1)+"/run-"+str(i+1)+"/rundata.csv")
        else:
            with open('rundata.csv','w',newline='')as file:
                writer = csv.writer(file)
                writer.writerow(['NO',prop,restart,inputs.cmprss,inputs.method,int(repeats/nodes),'NO'])
                writer.writerow(inputs.systems.values())
                writer.writerow(inputs.parameters.values())
                writer.writerow(inputs.Train.values())
                writer.writerow(inputs.clone.values())
                writer.writerow(inputs.paramz.values())
                writer.writerow(inham.EL.values())
                writer.writerow(inputs.prop.values())
                if(inputs.systems['System']=='MP'):
                    writer.writerow(inham.MP.values())
                elif(inputs.systems['System']=='HP'):
                    writer.writerow(inham.HP.values())
                else:
                    writer.writerow(inham.SB.values())
            for i in range (nodes):
                p=1
                q=1
                for file in glob.glob('run-'+str(i+1)+'/Outbs-00*.out'):
                    shutil.copy2(file,'../Outbsbackup/Outbs-'+str(p)+'.out_'+str(i+1)+'.out')
                    p=p+1
                if((inputs.clone['Cloning'])!='no'):
                    for file in glob.glob('run-'+str(i+1)+'/clonearr-*.out'):
                        shutil.copy2(file,'../Outbsbackup/clonearr-'+str(p)+'.out_'+str(i+1)+'.out')
                        q=q+1   
                shutil.copy2("rundata.csv","run-"+str(i+1))


    #If on a SGE machine make job submission file
    if(HPCFLG==1):
        number=random.randint(99999,1000000)
        file1="MCE"+str(number)+".sh"
        f=open(file1,"w")
        f.write("#$ -cwd -V \n")
        if(cores!=1):
            f.write("#$ -pe smp "+str(cores)+" \n") #Use shared memory parallel environemnt 
        f.write("#$ -l h_rt=40:00:00 \n")
        f.write("#$ -l h_vmem=4G \n")
        f.write("#$ -t 1-"+str(nodes)+" \n")
        f.write("date \n")
        f.write("cd "+EXDIR1+"/run-$SGE_TASK_ID/ \n")
        f.write("echo "'"Running on $HOSTNAME in folder $PWD" \n')
        f.write("module load mkl \n")
        f.write("time ./MCE.exe \n")
        f.write("date \n")
        f.close()
        if(cores!=1):
            os.environ["OMP_NUM_THREADS"]=str(cores)
        subprocess.call(['qsub',file1])

    else:
        if(cores!=1):
            os.environ["OMP_NUM_THREADS"]=str(cores)
        for i in range(nodes):
            SUBDIR=EXDIR1+"/run-"+str(i+1)
            subprocess.Popen('',executable=SUBDIR+"/MCE.exe",cwd=SUBDIR)
        
        

    