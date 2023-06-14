import numpy as np
from matplotlib import pyplot as plt
import os 
import shutil
import subprocess
import run
import glob
import sys
import inputs

time_end=inputs.prop['time_end']
time_start=inputs.prop['time_start']
nbf = inputs.nbf
repeats = run.repeats
nodes = run.nodes
repeats =int(repeats/nodes)

size=len(time_end)
time_end=float(time_end[:size-3])*10**1
dtinit=inputs.prop['dtinit']
size=len(dtinit)
dtinit=float(dtinit[:-3])*10**-3
time_start=0
incr= 0
no_t_step=int((time_end-time_start)/dtinit)



clone_freq=float(inputs.clone['clon_freq'])

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

if(inputs.clone["Cloning"]!='V1'):
    print("Only works with version 1 cloning")
else:
    print('version 1 condensing')







# for i in range(nodes):
#     SUBDIR=EXDIR1+"/run-"+str(i+1)
#     clonedir= EXDIR1+"/run-"+str(i+1)+"/clonetag.out"
#     x, y, z, a, b = np.loadtxt(clonedir, unpack=True)
#     normp=np.zeros((no_t_step+4,13))
#     # try:
#     timestep=z.astype(int)
#     SUBDIR=EXDIR1+"/run-"+str(i+1)

#     try:
#         for j in range(0,len(y)):
#             fp1=open(EXDIR1+"/run-"+str(i+1)+"/normpop-000"+str(int(y[j]))+".out")
#             for k, line in enumerate(fp1):
#                 if k >= (timestep[j]+4):
#                     list = line.split()
#                     list = [float(p) for p in list]
#                     list1=np.array(list)
#                     normp[k,1:]= np.add(normp[k,1:],list1[1:])    
#             fp1.close()

#     except TypeError:
#          print('here')
#          name= EXDIR1+"/run-"+str(i+1)+"/normpop-000"+str(int(y))+".out"
#          print(name)
#          fp2=open(EXDIR1+"/run-"+str(i+1)+"/normpop-000"+str(int(y))+".out")
#          for k, line in enumerate(fp2):
#             if k >= (timestep+4):
#                 list = line.split()
#                 list = [float(p) for p in list]
#                 list1=np.array(list)
#                 normp[k,1:]= np.add(normp[k,1:],list1[1:])    
#          fp2.close()
    

#     for j in range(0, repeats):
#         fp3=open(EXDIR1+"/run-"+str(i+1)+"/normpop-000"+str(int(j+1))+".out")
#         if(j==0):
#             for k, line in enumerate(fp3):
#                 if k>= 3:
#                     list = line.split()
#                     list = [float(p) for p in list]
#                     list1=np.array(list)
#                     normp[k,:]= np.add(normp[k,:],list1[:])
#         else:
#             for k, line in enumerate(fp3):
#                 if k>= 3:
#                     list = line.split()
#                     list = [float(p) for p in list]
#                     list1=np.array(list)
#                     normp[k,1:]= np.add(normp[k,1:],list1[1:])
#         fp3.close()
     
    
for i in range(nodes):
    SUBDIR=EXDIR1+"/run-"+str(i+1)
    clonedir= EXDIR1+"/run-"+str(i+1)+"/clonetag.out"
    parent, child, time_clone, normWP, normWC = np.loadtxt(clonedir, unpack=True)
    normp=np.zeros((no_t_step+4,13))
    timestep_clone = np.rint((time_clone/time_end)*no_t_step)
    totalreps= int(child[-1])
    shift= np.zeros((totalreps,3))
    for h in range(1,totalreps+1):
        # print(h)
        shift[h-1,0] = h
        if h in child:
            place = np.where(child==h)
            print('place is, ', place[0])
            place = int(place[0])
            shift[h-1,1] = timestep_clone[place]
            shift[h-1,2] = int(parent[place])  
    normweighting = np.ones((totalreps,no_t_step+1))
    # for h in range(1,repeats+1):
    #     for s in range(100,no_t_step+100,100):
    #         dir ='Outbs-'+str(h).zfill(4)+"-"+str(int(s)).zfill(4)
    #         # print(dir)
    # for h in range(0,totalreps-repeats):
    #     for s in range(round(int(timestep_clone[h]),-2),no_t_step-1,100):
    #         dir ='Outbs-'+str(int(child[h])).zfill(4)+"-"+str(int(s)+100-round(int(timestep_clone[h]),-2)).zfill(4)
    #         print(dir)
    for k in range(0,len(child)):
        row1 = int(parent[k]-1)
        row2 = int(child[k]-1)
        d = int(timestep_clone[k-1])
        hold = normweighting[row1, d]
        for m in range(int(timestep_clone[k]-1),no_t_step+1):
            normweighting[row1,m] = normweighting[row1, m] * normWP[k]
            normweighting[row2,m] = hold * normWC[k]
        for j in range(0,int(timestep_clone[k]-1)):
            normweighting[row2,j] = 0 
    for j in range(0,int(child[-1])):
        print(EXDIR1+"/run-"+str(i+1)+"/normpop-"+str(int(j+1)).zfill(4)+".out")
        fp1=open(EXDIR1+"/run-"+str(i+1)+"/normpop-"+str(int(j+1)).zfill(4)+".out")
        if(j==0):
            for k, line in enumerate(fp1):
                if k>= 3:
                    list = line.split()
                    list = [float(p) for p in list]
                    list1=np.array(list)
                    normp[k,0]= np.add(normp[k,0],list1[0])
                    list1 = list1 * normweighting[j,k-3]
                    normp[k,1:]= np.add(normp[k,1:],list1[1:])
        else:        
            for k, line in enumerate(fp1):
                if k >= 3:
                    list = line.split()
                    list = [float(p) for p in list]
                    list1=np.array(list)
                    list1= list1 *normweighting[j,k-3]
                    normp[k,1:]= np.add(normp[k,1:],list1[1:])    
        fp1.close()
    normp[:,1:]=normp[:,1:]/repeats

    # print(repeats)
    # print(normp[2003,:])
    
    # normp=normp[0,:]*repeats
    # print(normp[2003:])
    header = "Time Norm Re(ACF(t)) Im(ACF(t)) |ACF(t)| Re(Extra) Im(Extra) |Extra| Sum(HEhr) Pop1 Pop2 Pop1+Pop2 Pop2-Pop1 \n \n"
    np.savetxt(EXDIR1+"/run-"+str(i+1)+"/normpop.out",normp[3:,:], delimiter=' ', header=  header, comments= '')
    # plt.plot(normp[3:,0],normp[3:,12])
    # plt.show()