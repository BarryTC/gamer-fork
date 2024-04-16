#!/usr/bin/env python3.9

import sys
import math
import time
import os
import numpy as np
from datetime import datetime
from numpy import loadtxt
print(str(datetime.now()))


#-----------------------------------------------------------------------------------------------------------------
# user-specified parameters
fix_Delta_t           = float(0.2)           # data dump time interval [t_cross]
t_Start               = float(0.0)           # Input__DumpTable start time [t_cross]
t_End                 = float(10.0)          # Input__DumpTable end time [t_cross]

#-----------------------------------------------------------------------------------------------------------------
# return the numerical value (in str type) of "Parameter" in "FineName"
def Input_Parameter_Readline(FineName, Parameter, extraction_index):
   FileOpen = open('%s'%FineName, 'r').readlines()
   for Line_Now in FileOpen:
      if Parameter in Line_Now:
         return Line_Now.split()[extraction_index]
#Input_TestProb        = "../PAR_IC__Input__TestProb"
#t_cross               = float(Input_Parameter_Readline(Input_TestProb, "t_cross    ", 1))
NFW_concentration     = 10.0
t_cross               = NFW_concentration**1.5                      # the crossing time
t_cross_List          = np.linspace(t_Start, t_End, num=int(abs(t_End-t_Start)/fix_Delta_t)+1)
t_CodeUnit_List       = t_cross_List*t_cross

# output time table in [t_cross]
with open("Input__DumpTable_tcross", 'w') as master_output:
    # output header
    master_output.write("#Dump ID             Dump Time\n")
    # output simulation data
    for t_idx, t_Now in enumerate(t_cross_List):
            master_output.write('{:>8}'.format(t_idx))
            master_output.write('{:>22}'.format("%.2f"%t_Now))
            master_output.write("\n")
    master_output.write("***************END LINE***************\n")
print("Input__DumpTable_tcross COMPLETE")

# output time table in [Code_T]
with open("Input__DumpTable", 'w') as master_output:
    # output header
    master_output.write("#Dump ID             Dump Time\n")
    # output simulation data
    for t_idx, t_Now in enumerate(t_CodeUnit_List):
            master_output.write('{:>8}'.format(t_idx))
            master_output.write('{:>22}'.format("%.12e"%t_Now))
            master_output.write("\n")
    master_output.write("***************END LINE***************\n")
print("Input__DumpTable COMPLETE")
