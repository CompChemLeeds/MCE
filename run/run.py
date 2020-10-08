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
#   To run the program variables must be set/checked in input.py, inham.py and prop.py.
#   The following arguemtns then have ot be set in the run folder
#       1) The number of repeats
#       2) The number of folders/nodes
#       3) The number of parallel threads per folder/node (max8)
#
#
#
#
#
#
#########################################################################################
import sys
import socket
import prop
import inham
import input

#########################################################################################
#                              VARIABLES TO SET FOR SIMULATION                          #
#########################################################################################

# Number of repeats 
repeats=40
# Number of nodes/folders
nodes=1
#Number of parallel threads per folder/node (max 8)
threads=4
# Name of running folder 
# Default : <method>-<system>-<random number> ie CCS-HP-31254
# Otherwise:  <method>-<system>-<runfolder string>
Runfolder='testforMWE'
# Generate Basis Set? YES/NO
gen='YES'
# Propagate Basis Set? YES/NO
prop='YES'
# Restart? YES/NO
# To restart a timedout run set to yes and rerun this script from the execution folder
restart='NO'
#########################################################################################
#                                   END OF INPUTS                                       #
#########################################################################################
#                * NO NEED TO SCROLL FURTHER IF USING AS BLACKBOX *                     #
#########################################################################################

if(isinstance(repeats,int)==False):
    sys.exit("Number of repeats must be an integer")
elif(isinstance(nodes,int)==False):
    sys.exit("Number of folders must be an integer")
elif(isinstance(threads,int)==False):
    sys.exit("Number of parallel threads must be an integer")
elif(repeats<1):
    sys.exit("Not enough runs selected. Must be 1 or greater")
elif(nodes<1):
    sys.exit("Not enough nodes selected. Must be 1 or greater")
elif(nodes>100):
    sys.exit("Too many nodes. Maximum of 100 simultaneous submisions")
elif(threads>8):
    sys.exit("Too many threads selected. Maximum of 8 available")
elif(threads<1):
    sys.exit("Not enough threads selected. Must be 1 or greater")
elif((repeats/nodes)>5000):
    sys.exit("Too many repeats per folder. Must be less than 5000")
elif((input.Conjugate_Repeats=='YES')and((repeats%(2*nodes*threads))!=0)):
    sys.exit("Number of repeats not valid for conjugate repetition. Should be integer multiple of 2*threads*nodes")
elif((repeats%(nodes*threads))!=0):
    sys.exit("Number of repeats must be an integer multiple of cores*folders")
elif(nodes*threads>100):
    sys.exit("Total number of cores should stay below 100")
else:
    print("Arguments checked")    
    if(socket.gethostname()==("arc3"or"arc4")):
        HPCFLG=1
    else:
        HPCFLG=0

print(HPCFLG)