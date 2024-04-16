import argparse
import sys
import yt
import linecache
import numpy as np
import math
import os
import matplotlib.image
import matplotlib.pyplot as plt
from numpy import loadtxt


#-------------------------------------------------------------------------------------------------------------------------
# user-specified parameters
unit_mode             = 1      # use M_h/R_scale/G [1] or astrophysical units Msun/kpc/Gyr [2]
vmax_fac              = 1.5    # prefactor that adjusts the value of v_max in the 3D velocity dispersion plot; default is set to [1]
dpi                   = 110


#-------------------------------------------------------------------------------------------------------------------------
# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the gas slices and projections' )

parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix',
                     help='path prefix [%(default)s]', default='../' )
parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )

args=parser.parse_args()

idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix


#-------------------------------------------------------------------------------------------------------------------------
# return the numerical value (in str type) of "Parameter" in "FineName"
def Input_Parameter_Readline(FineName, Parameter, extraction_index):
   FileOpen = open('%s'%FineName, 'r').readlines()
   for Line_Now in FileOpen:
      if Parameter in Line_Now:
         return Line_Now.split()[extraction_index]
# find the nearest value
# https://stackoverflow.com/a/2566508
def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

# load parameters
Input_TestProb        = "../../PAR_IC__Input__TestProb"
NFW_concentration     = float(Input_Parameter_Readline(Input_TestProb, "NFW_concentration    ", 1))
BoxSize_Half          = float(Input_Parameter_Readline(Input_TestProb, "BOX_SIZE    ", 1))/2.0
Radial_Bin_End        = math.log10(BoxSize_Half)
Radial_Bin_Start      = Radial_Bin_End-3.8

# load analytical/expected anisotropy table
anisotropy_filename   = "../Radial_Anisotropy_Analysis/Particle_Radial_Anisotropy_000000" # input master list filename
anisotropy_data       = loadtxt("%s"%anisotropy_filename, skiprows=1, dtype=float)
analytical_r          = anisotropy_data[:,0].astype(float)
analytical_v3D        = anisotropy_data[:,2].astype(float)
analytical_rho        = anisotropy_data[:,3].astype(float)
analytical_mass       = anisotropy_data[:,5].astype(float)
analytical_par_count_per_bin = anisotropy_data[:,4].astype(int)
# remove bins with "N_par < 10"
analytical_r          = analytical_r[analytical_par_count_per_bin >= 10]
analytical_v3D        = analytical_v3D[analytical_par_count_per_bin >= 10]
analytical_rho        = analytical_rho[analytical_par_count_per_bin >= 10]
analytical_mass       = analytical_mass[analytical_par_count_per_bin >= 10]
analytical_beta       = analytical_r*0.0 # isotropic

# set the plot axis limits based on the initial profile
xlim_lower            = 0.85*analytical_r[0]
xlim_upper            = 2.2*analytical_r[-1]
vlim_upper            = 1.5*np.max(analytical_v3D)
rholim_upper          = 5.8*np.max(analytical_v3D)
rholim_lower          = 1.0e-13*np.max(analytical_v3D)

# load CM position and velocity
CM_PosVel_data        = loadtxt("../../Record__Center", skiprows=1, dtype=float)
CM_Time_list          = CM_PosVel_data[:,0].astype(float)
NPar_Bound_list       = CM_PosVel_data[:,2].astype(int)
CM_x_list             = CM_PosVel_data[:,4].astype(float)
CM_y_list             = CM_PosVel_data[:,5].astype(float)
CM_z_list             = CM_PosVel_data[:,6].astype(float)
CM_Vx_list            = CM_PosVel_data[:,7].astype(float)
CM_Vy_list            = CM_PosVel_data[:,8].astype(float)
CM_Vz_list            = CM_PosVel_data[:,9].astype(float)

# import "Input__DumpTable"
Input__DumpTable_data = loadtxt("../../Input__DumpTable", skiprows=1, dtype=str)
Data_ID_List          = Input__DumpTable_data[:-1,0].astype(int)
Time_List             = Input__DumpTable_data[:-1,1].astype(float)
Input__DumpTable_data = loadtxt("../../Input__DumpTable_tcross", skiprows=1, dtype=str)
Time_List_tcross      = Input__DumpTable_data[:-1,1].astype(float)

#-------------------------------------------------------------------------------------------------------------------------

for idx in range(idx_start,idx_end+1,didx):
   print("Processing idx = %i"%idx)
   # subhalo CM position
   Data_idx           = np.where(Data_ID_List == idx)[0][0]
   Data_idx_2nd       = np.nonzero(np.isclose(Data_ID_List, idx))[0][0]
   Time_Now           = Time_List[Data_idx_2nd]
   Time_Now_tcross    = Time_List_tcross[Data_idx_2nd]
   if (Data_idx != Data_idx_2nd):
      print("Check1: Frame_idx (%i) != np.isclose (%i)"%(Frame_idx, Frame_idx_2nd))
   Frame_idx          = np.where(CM_Time_list == Time_Now)[0][0]
   Frame_idx_2nd      = np.nonzero(np.isclose(CM_Time_list, Time_Now))[0][0]
   Coord_CoM_Pos_Vec  = np.array([CM_x_list[Frame_idx], CM_y_list[Frame_idx], CM_z_list[Frame_idx]])-BoxSize_Half
   Vel_CoM_Pos_Vec    = np.array([CM_Vx_list[Frame_idx], CM_Vy_list[Frame_idx], CM_Vz_list[Frame_idx]])
   if (Frame_idx != Frame_idx_2nd):
      print("Check2: Frame_idx (%i) != np.isclose (%i)"%(Frame_idx, Frame_idx_2nd))

   #---------------------------------------------------------------------------------------------------------------------
   #---------------------------------------------------------------------------------------------------------------------
   # load GAMER anisotropy table
   anisotropy_filename    = "../Radial_Anisotropy_Analysis/Particle_Radial_Anisotropy_%06d"%idx # input master list filename
   data                   = loadtxt("%s"%anisotropy_filename, skiprows=1, dtype=float)
   processed_r            = data[:,0].astype(float)
   processed_beta         = data[:,1].astype(float)
   processed_v3D          = data[:,2].astype(float)
   par_count_per_bin      = data[:,4].astype(int)
   processed_mass         = data[:,5].astype(float)
   processed_rho          = data[:,3].astype(float)
   # remove bins with "N_par < 10"
   processed_r            = processed_r[par_count_per_bin >= 10]
   processed_beta         = processed_beta[par_count_per_bin >= 10]
   processed_v3D          = processed_v3D[par_count_per_bin >= 10]
   processed_mass         = processed_mass[par_count_per_bin >= 10]
   processed_rho          = processed_rho[par_count_per_bin >= 10]
   par_count_per_bin      = par_count_per_bin[par_count_per_bin >= 10]

   # plot shell-averaged physical profiles
   # plot the density profile
   fig, axs = plt.subplots(5, figsize=(6, 13))
   axs[0].loglog(analytical_r, analytical_rho, linestyle='-.', color='orange', label = "PAR_IC", zorder=1)
   axs[0].loglog(processed_r, processed_rho, alpha=0.175, zorder=5)
   axs[0].scatter(processed_r, processed_rho, s=12, facecolors='none', edgecolors='C0', label = "GAMER", zorder=10)
   axs[0].set_xlim(xlim_lower, xlim_upper)
   axs[0].set_ylim(rholim_lower, rholim_upper)
   axs[0].set_xlabel("r/$r_\mathrm{s}$", fontsize=12)
   axs[0].set_ylabel(r"$\rho [m_0/r_\mathrm{s}^3]$", fontsize=12)
   axs[0].legend(loc="lower left", fontsize=12)
   axs[0].set_title("Radial Density Profile at $t = %.2f t_\mathrm{cross}$"%Time_Now_tcross, fontsize=12)
   axs[0].tick_params(axis='x', which='major', direction='in')
   axs[0].tick_params(axis='x', which='minor', direction='in')
   axs[0].tick_params(axis='y', which='major', direction='in')
   axs[0].tick_params(axis='y', which='minor', direction='in')
   axs[0].xaxis.set_ticks_position('both')
   axs[0].yaxis.set_ticks_position('both')
   # plot the enclosed mass profile
   axs[1].loglog(analytical_r, analytical_mass, linestyle='-.', color='orange', label = "PAR_IC", zorder=1)
   axs[1].loglog(processed_r, processed_mass, alpha=0.175, zorder=5)
   axs[1].scatter(processed_r, processed_mass, s=12, facecolors='none', edgecolors='C0', label = "GAMER", zorder=10)
   axs[1].set_xlim(xlim_lower, xlim_upper)
   axs[1].set_ylim(1.5e-5,2.8)
   axs[1].set_xlabel("r/$r_\mathrm{s}$", fontsize=12)
   axs[1].set_ylabel("Enclosed Mass [$m_0$]", fontsize=12)
   axs[1].legend(loc="lower right", fontsize=12)
   axs[1].set_title("Enclosed Mass Profile", fontsize=12)
   axs[1].tick_params(axis='x', which='major', direction='in')
   axs[1].tick_params(axis='x', which='minor', direction='in')
   axs[1].tick_params(axis='y', which='major', direction='in')
   axs[1].tick_params(axis='y', which='minor', direction='in')
   axs[1].xaxis.set_ticks_position('both')
   axs[1].yaxis.set_ticks_position('both')
   # plot the radial anisotropy profile
   axs[2].semilogx(analytical_r, analytical_beta, linestyle='-.', color='orange', label = "Isotropic", zorder=0)
   axs[2].semilogx(processed_r, processed_beta, alpha=0.175, zorder=5)
   axs[2].scatter(processed_r, processed_beta, s=12, facecolors='none', edgecolors='C0', label = "GAMER", zorder=10)
   axs[2].set_xlim(xlim_lower, xlim_upper)
   axs[2].set_ylim(-0.7, 1.25)
   axs[2].set_xlabel("r/$r_\mathrm{s}$", fontsize=12)
   axs[2].set_ylabel(r"Radial Anisotropy $\beta$", fontsize=12)
   axs[2].set_title("Radial Anisotropy (Isotropic)", fontsize=12)
   axs[2].legend(loc="upper left", fontsize=12)
   axs[2].tick_params(axis='x', which='major', direction='in')
   axs[2].tick_params(axis='x', which='minor', direction='in')
   axs[2].tick_params(axis='y', which='major', direction='in')
   axs[2].tick_params(axis='y', which='minor', direction='in')
   axs[2].xaxis.set_ticks_position('both')
   axs[2].yaxis.set_ticks_position('both')
   # plot the 3D velocity profile
   axs[3].semilogx(analytical_r, analytical_v3D, linestyle='-.', color='orange', label = "PAR_IC", zorder=0)
   axs[3].semilogx(processed_r, processed_v3D, alpha=0.175, zorder=5)
   axs[3].scatter(processed_r, processed_v3D, s=12, facecolors='none', edgecolors='C0', label = "GAMER", zorder=10)
   axs[3].set_xlim(xlim_lower, xlim_upper)
   axs[3].set_ylim(0,vlim_upper)
   axs[3].set_xlabel("r/$r_\mathrm{s}$", fontsize=12)
   axs[3].set_ylabel("$\sigma_\mathrm{3D}$ [$G = 1$]", fontsize=12)
   axs[3].legend(loc="upper right", fontsize=12)
   axs[3].set_title("3D Velocity Dispersion", fontsize=12)
   axs[3].tick_params(axis='x', which='major', direction='in')
   axs[3].tick_params(axis='x', which='minor', direction='in')
   axs[3].tick_params(axis='y', which='major', direction='in')
   axs[3].tick_params(axis='y', which='minor', direction='in')
   axs[3].xaxis.set_ticks_position('both')
   axs[3].yaxis.set_ticks_position('both')
   # plot the total particle counts in each radial bin
   axs[4].loglog(processed_r, par_count_per_bin, alpha=0.175, zorder=5)
   axs[4].scatter(processed_r, par_count_per_bin, s=12, facecolors='none', edgecolors='C0', label = "GAMER", zorder=10)
   axs[4].set_xlim(xlim_lower, xlim_upper)
   axs[4].set_ylim(0.85, 2.5e6)
   axs[4].set_xlabel("r/$r_\mathrm{s}$", fontsize=12)
   axs[4].set_ylabel("$N_\mathrm{particle}$", fontsize=12)
   axs[4].legend(loc="upper right", fontsize=12)
   axs[4].set_title("Particle Counts Per Bin", fontsize=12)
   axs[4].tick_params(axis='x', which='major', direction='in')
   axs[4].tick_params(axis='x', which='minor', direction='in')
   axs[4].tick_params(axis='y', which='major', direction='in')
   axs[4].tick_params(axis='y', which='minor', direction='in')
   axs[4].xaxis.set_ticks_position('both')
   axs[4].yaxis.set_ticks_position('both')
   fig.tight_layout()
   plt.savefig("Data_%06d_Density_Anisotropy_Profiles.png"%idx, dpi = dpi)
   plt.close()
