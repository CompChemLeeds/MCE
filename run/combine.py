# Data combination script 
# Written by O. A. Bramley 04/01/21
#
# This is based on a bash cript written by C. Symonds. It combines the normpop files from
# multiple partial runs. The script calls an averaging program subavrg.exe for static 
# stepsizes and the interpolation program, interpolate.ex to produce a final normpop.out
# file for each folder as would be needed by the collate.py script which is called at the 
# end.

import run

run.nodes