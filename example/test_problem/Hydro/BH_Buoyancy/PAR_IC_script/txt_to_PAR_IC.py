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
NFW_concentration      = 10.0
Input_PAR_IC_Filename  = "./snapshot.txt"
Boxsize                = 3000.0                             # [NFW_Rscale]
Center_of_Mass_X       = 0.5*Boxsize                        # central x coordinate, i.e. simulation box center [NFW_Rscale]
Center_of_Mass_Y       = 0.5*Boxsize                        # central y coordinate, i.e. simulation box center [NFW_Rscale]
Center_of_Mass_Z       = 0.5*Boxsize                        # central z coordinate, i.e. simulation box center [NFW_Rscale]
Par_type               = 1.0


#-----------------------------------------------------------------------------------------------------------------
# constants and conversions in model units "G = NFW_Rscale = NFW_Mvir = 1"; see "van den Bosch1 and Ogiya, MNRAS (2018)"
U_G           = 1.0                                         # gravitational constant
NFW_Rscale    = 1.0
NFW_Mvir      = 1.0
fNFW          = math.log(1.0+NFW_concentration)-NFW_concentration/(1.0+NFW_concentration)
NFW_rho0      = NFW_Mvir/(4.0*math.pi*NFW_Rscale**3.0*fNFW)
NFW_Vvir      = 1.0/(NFW_concentration)**0.5                # the initial virial velocity
t_cross       = NFW_concentration**1.5                      # the crossing time
U_Gyr         = t_cross/2.006                               # gigayear
U_yr          = U_Gyr/1.0e9                                 # tropical year
U_Myr         = U_yr/1.0e6                                  # megayear
U_sec         = U_yr/31556925.2
# by fixing "NFW_Rscale = 1.0*U_cm" for simplicity
U_cm          = 1.0*NFW_Rscale
U_gram        = 6.6738*1e-8*U_cm**3/(U_G*U_sec**2)
U_pc          = U_cm*3.08567758149*1e18                     # parsec
U_kpc         = 1.0e3*U_pc
U_Mpc         = 1.0e6*U_pc
U_Msun        = U_gram*1.9885*1e33                          # solar mass
U_kg          = U_Msun/(1.9885*1e30)
U_meter       = U_pc /(3.08567758149*1e16)
U_km          = 1.0e3*U_meter
U_c           = 299792458*U_meter/U_sec                     # speed of light
U_eV          = 1.782661907*1e-36*U_kg*U_c**2               # electronvolt
U_MeV         = 1e6*U_eV

print("The current unit settings are: NFW_concentration = %.2f, UNIT_L/cm = %.10e, UNIT_M/gram = %.10e, UNIT_T/sec = %.10e"%(NFW_concentration, 1.0/U_cm, 1.0/U_gram, 1.0/U_sec))


#-----------------------------------------------------------------------------------------------------------------
# load input data PAR_IC text file
Input_PAR_IC_data     = loadtxt(Input_PAR_IC_Filename, skiprows=0, dtype=float)
Output_Mass           = (Input_PAR_IC_data[:,0].astype(float))
Output_Pos_x          = (Input_PAR_IC_data[:,1].astype(float)) + Center_of_Mass_X
Output_Pos_y          = (Input_PAR_IC_data[:,2].astype(float)) + Center_of_Mass_Y
Output_Pos_z          = (Input_PAR_IC_data[:,3].astype(float)) + Center_of_Mass_Z
Output_Vel_x          = (Input_PAR_IC_data[:,4].astype(float))
Output_Vel_y          = (Input_PAR_IC_data[:,5].astype(float))
Output_Vel_z          = (Input_PAR_IC_data[:,6].astype(float))
Output_ParType        = Output_Mass*0.0 + Par_type
Output_ParTag         = Output_Mass*0.0 + 1.0
Output_ParEnergy      = Output_Mass*0.0
Output_ParAngMom      = Output_Mass*0.0
Output_ParEnergyBac   = Output_Mass*0.0
NPar_AllRank          = len(Output_Mass)                    # total number of particles
Output_ParID          = np.array(range(NPar_AllRank))
print([Output_Mass[0], Output_Pos_x[0], Output_Vel_x[0]])
print("Maximum mass = %.14e [Unit_M]; Minimum mass = %.14e [Unit_M]"%(np.max(Output_Mass),np.min(Output_Mass)))
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
   writer.write("t_cross                 %.14e         # [UNIT_T]\n"%t_cross)
   writer.write("BOX_SIZE                %.14e         # simulation box size [UNIT_L]\n"%Boxsize)
   writer.write("NFW_concentration       %.14e         # NFW concentration parameter\n"%NFW_concentration)
   writer.write("NPar_AllRank            %d                      # total number of particles\n"%NPar_AllRank)
   writer.write("Par_type                %.14e         # GAMER input particle types (for particle tagging)\n"%Par_type)
   writer.write("Particle_Mass           %.14e         # individual particle mass [UNIT_M]\n"%(np.min(Output_Mass)))
   writer.write("Ini_R_vir               %.14e         # initial NFW virial radius [UNIT_L]\n"%(NFW_Rscale*NFW_concentration))
   writer.write("Ini_vel_cir             %.14e         # initial circular velocity at virial radius [UNIT_L/UNIT_T]\n"%NFW_Vvir)
   writer.write("Center_of_Mass_X        %.14e         # central x coordinate, i.e. simulation box center [UNIT_L]\n"%Center_of_Mass_X)
   writer.write("Center_of_Mass_Y        %.14e         # central x coordinate, i.e. simulation box center [UNIT_L]\n"%Center_of_Mass_Y)
   writer.write("Center_of_Mass_Z        %.14e         # central x coordinate, i.e. simulation box center [UNIT_L]\n"%Center_of_Mass_Z)
print('PAR_IC__Input__TestProb complete')
