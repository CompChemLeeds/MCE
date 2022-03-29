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
repeats = run.repeats
nodes = run.nodes
repeats =int(repeats/nodes)

size=len(time_end)
time_end=float(time_end[:size-3])*10**1
dtinit=inputs.prop['dtinit']
size=len(dtinit)
dtinit=float(dtinit[:-3])*10**-3
time_start=0

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

if(inputs.clone["Cloning"]!='v1'):
    print("Only works with version 1 cloning")
else:
    print('version 1 condensing')

# for i in range(nodes):
#     SUBDIR=EXDIR1+"/run-"+str(i+1)
#     clonedir= EXDIR1+"/run-"+str(i+1)+"/clonetag.out"
#     x, y, z, a, b = np.loadtxt(clonedir, unpack=True)
#     x = list(x)
#     y = list(y) 
#     z = list(z)
#     a = list(a)
#     b = list(b)
#     print(x,y,z,a,b)
#     totalreps= int(y[-1])
#     normweighting = np.ones((totalreps,no_t_step))
#     for k in range(0,len(y)):
#         row1 = int(x[k]-1)
#         row2 = int(y[k]-1)
#         d = int(z[k-1])
#         hold = normweighting[row1, d]
#         for m in range(int(z[k]-1),no_t_step):
#             normweighting[row1,m] = normweighting[row1, m] * a[k]
#             normweighting[row2,m] = hold * b[k]
#         for j in range(0,int(z[k]-1)):
#             normweighting[row2,j] = 0 
# print(normweighting[0:,0:10])
# print(normweighting[0:,1490:1510])
# print(normweighting[0:,1990:])






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
    x, y, z, a, b = np.loadtxt(clonedir, unpack=True)
    normp=np.zeros((no_t_step+4,13))
    x = list(x)
    y = list(y) 
    z = list(z)
    a = list(a)
    b = list(b)
    print(x,y,z,a,b)
    totalreps= int(y[-1])
    normweighting = np.ones((totalreps,no_t_step+1))
    for k in range(0,len(y)):
        row1 = int(x[k]-1)
        row2 = int(y[k]-1)
        d = int(z[k-1])
        hold = normweighting[row1, d]
        for m in range(int(z[k]-1),no_t_step+1):
            normweighting[row1,m] = normweighting[row1, m] * a[k]
            normweighting[row2,m] = hold * b[k]
        for j in range(0,int(z[k]-1)):
            normweighting[row2,j] = 0 
    # timestep=z.astype(int)
    SUBDIR=EXDIR1+"/run-"+str(i+1)
    for j in range(0,int(y[-1])):
        print(EXDIR1+"/run-"+str(i+1)+"/normpop-"+str(int(j+1)).zfill(4)+".out")
        fp1=open(EXDIR1+"/run-"+str(i+1)+"/normpop-"+str(int(j+1)).zfill(4)+".out")
        if(j==0):
            for k, line in enumerate(fp1):
                if k>= 3:
                    list = line.split()
                    list = [float(p) for p in list]
                    list1=np.array(list)
                    # print(k, list1)
                    normp[k,0]= np.add(normp[k,0],list1[0])
                    list1 = list1 * normweighting[j,k-3]
                    # print(k, list1)
                    normp[k,1:]= np.add(normp[k,1:],list1[1:])
        else:        
            for k, line in enumerate(fp1):
                if k >= 3:
                    list = line.split()
                    list = [float(p) for p in list]
                    # print(EXDIR1+"/run-"+str(i+1)+"/normpop-000"+str(int(j+1))+".out")
                    # print(k-3, normweighting[j,0:10])
                    list1=np.array(list)
                    list1= list1 *normweighting[j,k-3]
                    normp[k,1:]= np.add(normp[k,1:],list1[1:])    
        fp1.close()


    # print(repeats)
    # print(normp[2003,:])
    normp[:,1:]=normp[:,1:]/repeats
    # normp=normp[0,:]*repeats
    # print(normp[2003:])
    header = "Time Norm Re(ACF(t)) Im(ACF(t)) |ACF(t)| Re(Extra) Im(Extra) |Extra| Sum(HEhr) Pop1 Pop2 Pop1+Pop2 Pop2-Pop1 \n \n"
    np.savetxt(EXDIR1+"/run-"+str(i+1)+"/normpop.out",normp[3:,:], delimiter=' ', header=  header, comments= '')
    # plt.plot(normp[3:,0],normp[3:,12])
    # plt.show()