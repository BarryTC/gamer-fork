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
Inputfile_DM           = "./DM_IC_List"
Inputfile_Stars        = "./Stars_IC_List"
Boxsize                = 200.0                              # [kpc]
Center_of_Mass_X       = 0.5*Boxsize                        # central x coordinate, i.e. simulation box center [NFW_Rscale]
Center_of_Mass_Y       = 0.5*Boxsize                        # central y coordinate, i.e. simulation box center [NFW_Rscale]
Center_of_Mass_Z       = 0.5*Boxsize                        # central z coordinate, i.e. simulation box center [NFW_Rscale]
DM_ParType             = 2.0                                # https://github.com/gamer-project/gamer/wiki/Initial-Conditions#particles-1
Stars_ParType          = 3.0


#-----------------------------------------------------------------------------------------------------------------
# constants and conversions in astrophysical units {Msun, pc, Myr} = {1, 1, 1};
U_Msun   = 1.0                                   # solar mass
U_pc     = 1.0                                   # parsec
U_Myr    = 1.0                                   # megayear
U_kpc    = 1.0e3*U_pc
U_Mpc    = 1.0e6*U_pc
U_gram   = U_Msun/(1.9885*1e33)
U_kg     = U_Msun/(1.9885*1e30)
U_cm     = U_pc /(3.08567758149*1e18)
U_meter  = U_pc /(3.08567758149*1e16)
U_km     = 1.0e3*U_meter
U_yr     = U_Myr/1.0e6
U_Gyr    = 1.0e9*U_yr
U_sec    = U_yr/31556925.2                       # tropical year
U_c      = 299792458*U_meter/U_sec               # speed of light
U_eV     = 1.782661907*1e-36*U_kg*U_c**2         # electronvolt
U_MeV    = 1.0e6*U_eV
U_G      = 6.6738*1e-8*U_cm**3/(U_gram*U_sec**2) # gravitational constant
U_hbar   = 6.582119514*1e-22*U_MeV*U_sec
U_V_conv = (U_pc/U_Myr)/(U_km/U_sec)


#-----------------------------------------------------------------------------------------------------------------
# load input data PAR_IC text file
Input_PAR_IC_DM       = loadtxt(Inputfile_DM, skiprows=1, dtype=float)
Output_Mass_DM        = (Input_PAR_IC_DM[:,0].astype(float))
Output_Pos_x_DM       = ((Input_PAR_IC_DM[:,1].astype(float)) + Center_of_Mass_X)*U_kpc
Output_Pos_y_DM       = ((Input_PAR_IC_DM[:,2].astype(float)) + Center_of_Mass_Y)*U_kpc
Output_Pos_z_DM       = ((Input_PAR_IC_DM[:,3].astype(float)) + Center_of_Mass_Z)*U_kpc
Output_Vel_x_DM       = (Input_PAR_IC_DM[:,4].astype(float))*U_V_conv
Output_Vel_y_DM       = (Input_PAR_IC_DM[:,5].astype(float))*U_V_conv
Output_Vel_z_DM       = (Input_PAR_IC_DM[:,6].astype(float))*U_V_conv
Output_ParType_DM     = Output_Mass_DM*0.0 + DM_ParType

Input_PAR_IC_Stars    = loadtxt(Inputfile_Stars, skiprows=1, dtype=float)
Output_Mass_Stars     = (Input_PAR_IC_Stars[:,0].astype(float))
Output_Pos_x_Stars    = ((Input_PAR_IC_Stars[:,1].astype(float)) + Center_of_Mass_X)*U_kpc
Output_Pos_y_Stars    = ((Input_PAR_IC_Stars[:,2].astype(float)) + Center_of_Mass_Y)*U_kpc
Output_Pos_z_Stars    = ((Input_PAR_IC_Stars[:,3].astype(float)) + Center_of_Mass_Z)*U_kpc
Output_Vel_x_Stars    = (Input_PAR_IC_Stars[:,4].astype(float))*U_V_conv
Output_Vel_y_Stars    = (Input_PAR_IC_Stars[:,5].astype(float))*U_V_conv
Output_Vel_z_Stars    = (Input_PAR_IC_Stars[:,6].astype(float))*U_V_conv
Output_ParType_Stars  = Output_Mass_Stars*0.0 + Stars_ParType

NPar_AllRank          = len(Output_Mass_DM) + len(Output_Mass_Stars)     # total number of particles

Output_Mass           = np.concatenate((Output_Mass_DM, Output_Mass_Stars), axis=0)
Output_Pos_x          = np.concatenate((Output_Pos_x_DM, Output_Pos_x_Stars), axis=0)
Output_Pos_y          = np.concatenate((Output_Pos_y_DM, Output_Pos_y_Stars), axis=0)
Output_Pos_z          = np.concatenate((Output_Pos_z_DM, Output_Pos_z_Stars), axis=0)
Output_Vel_x          = np.concatenate((Output_Vel_x_DM, Output_Vel_x_Stars), axis=0)
Output_Vel_y          = np.concatenate((Output_Vel_y_DM, Output_Vel_y_Stars), axis=0)
Output_Vel_z          = np.concatenate((Output_Vel_z_DM, Output_Vel_z_Stars), axis=0)
Output_ParType        = np.concatenate((Output_ParType_DM, Output_ParType_Stars), axis=0)
Output_ParTag         = Output_Mass*0.0 + 1.0
Output_ParEnergy      = Output_Mass*0.0
Output_ParAngMom      = Output_Mass*0.0
Output_ParEnergyBac   = Output_Mass*0.0
Output_ParID          = np.array(range(NPar_AllRank))+1
print("Output a total number of %i particles"%NPar_AllRank)

# output final PAR_IC binary file
with open('PAR_IC', 'wb') as output:
       output.write(Output_Mass.astype('f').tobytes())
       output.write(Output_Pos_x.astype('f').tobytes())
       output.write(Output_Pos_y.astype('f').tobytes())
       output.write(Output_Pos_z.astype('f').tobytes())
       output.write(Output_Vel_x.astype('f').tobytes())
       output.write(Output_Vel_y.astype('f').tobytes())
       output.write(Output_Vel_z.astype('f').tobytes())
       output.write(Output_ParType.astype('f').tobytes())
       output.write(Output_ParTag.astype('f').tobytes())
       output.write(Output_ParEnergy.astype('f').tobytes())
       output.write(Output_ParAngMom.astype('f').tobytes())
       output.write(Output_ParEnergyBac.astype('f').tobytes())
       output.write(Output_ParID.astype('f').tobytes())
output.close()
print('PAR_IC complete')

# output adopted simulation parameters
with open('PAR_IC__Input__TestProb', 'w') as writer:
   writer.write("UNIT_L                  %.14e\n"%(1.0/U_cm))
   writer.write("UNIT_M                  %.14e\n"%(1.0/U_gram))
   writer.write("UNIT_T                  %.14e\n"%(1.0/U_sec))
   writer.write("BOX_SIZE                %.14e         # simulation box size [UNIT_L]\n"%Boxsize)
   writer.write("NPar_AllRank            %d            # total number of particles\n"%NPar_AllRank)
   writer.write("Center_of_Mass_X        %.14e         # central x coordinate, i.e. simulation box center [UNIT_L]\n"%(Center_of_Mass_X*U_kpc))
   writer.write("Center_of_Mass_Y        %.14e         # central x coordinate, i.e. simulation box center [UNIT_L]\n"%(Center_of_Mass_Y*U_kpc))
   writer.write("Center_of_Mass_Z        %.14e         # central x coordinate, i.e. simulation box center [UNIT_L]\n"%(Center_of_Mass_Z*U_kpc))
print('PAR_IC__Input__TestProb complete')
