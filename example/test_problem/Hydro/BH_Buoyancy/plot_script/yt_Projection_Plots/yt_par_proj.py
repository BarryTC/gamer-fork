import argparse
import sys
import yt
import numpy as np
import math
import os
import os.path
import matplotlib.image
import matplotlib.pyplot as plt
from numpy import loadtxt


#-------------------------------------------------------------------------------------------------------------------------
# user-specified parameters
vmax_fac        = 1.5    # prefactor that adjusts the value of v_max in the 3D velocity dispersion plot; default is set to [1]
mass_lim_upper  = 2.5e-2
mass_lim_lower  = 0.5e-6
zoom_factor     = 50
dpi             = 150

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
# NFW density profile
def NFW_Profile(rnow, rho_scale, r_scale):
   return rho_scale/((float(rnow)/r_scale)*(1.0+float(rnow)/r_scale)**2.0)
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
Input_TestProb    = "../../PAR_IC__Input__TestProb"
UNIT_L            = float(Input_Parameter_Readline(Input_TestProb, "UNIT_L    ", 1))
UNIT_M            = float(Input_Parameter_Readline(Input_TestProb, "UNIT_M    ", 1))
UNIT_T            = float(Input_Parameter_Readline(Input_TestProb, "UNIT_T    ", 1))
BOX_SIZE          = float(Input_Parameter_Readline(Input_TestProb, "BOX_SIZE    ", 1))
NFW_concentration = float(Input_Parameter_Readline(Input_TestProb, "NFW_concentration    ", 1))
NPar_AllRank      = float(Input_Parameter_Readline(Input_TestProb, "NPar_AllRank    ", 1))
Par_type          = float(Input_Parameter_Readline(Input_TestProb, "Par_type    ", 1))
Particle_Mass     = float(Input_Parameter_Readline(Input_TestProb, "Particle_Mass    ", 1))
Ini_R_vir         = float(Input_Parameter_Readline(Input_TestProb, "Ini_R_vir    ", 1))
Ini_vel_cir       = float(Input_Parameter_Readline(Input_TestProb, "Ini_vel_cir    ", 1))
Center_of_Mass_X  = float(Input_Parameter_Readline(Input_TestProb, "Center_of_Mass_X    ", 1))
Center_of_Mass_Y  = float(Input_Parameter_Readline(Input_TestProb, "Center_of_Mass_Y    ", 1))
Center_of_Mass_Z  = float(Input_Parameter_Readline(Input_TestProb, "Center_of_Mass_Z    ", 1))

# load CM position and velocity
CM_PosVel_data       = loadtxt("../../Record__Center", skiprows=1, dtype=float)
CM_Time_list         = CM_PosVel_data[:,0].astype(float)
NPar_Bound           = CM_PosVel_data[:,2].astype(float)
CM_x_list            = CM_PosVel_data[:,4].astype(float)
CM_y_list            = CM_PosVel_data[:,5].astype(float)
CM_z_list            = CM_PosVel_data[:,6].astype(float)
CM_Vx_list           = CM_PosVel_data[:,7].astype(float)
CM_Vy_list           = CM_PosVel_data[:,8].astype(float)
CM_Vz_list           = CM_PosVel_data[:,9].astype(float)
NPar_Bound_Max       = float(NPar_Bound.max())

# import "Input__DumpTable"
Input__DumpTable_data = loadtxt("../../Input__DumpTable", skiprows=1, dtype=str)
Data_ID_List          = Input__DumpTable_data[:-1,0].astype(int)
Time_List             = Input__DumpTable_data[:-1,1].astype(float)
Input__DumpTable_data = loadtxt("../../Input__DumpTable_tcross", skiprows=1, dtype=str)
Time_List_tcross      = Input__DumpTable_data[:-1,1].astype(float)


#-------------------------------------------------------------------------------------------------------------------------
field          = 'particle_mass'
for idx in range(idx_start,idx_end+1,didx):
   # load simulation data with yt
   ds          = yt.load('../../Data_%06d'%idx)
   Data_ID_Now = ds.parameters["DumpID"]
   BoxSize     = float(ds.domain_width[1].in_units("code_length"))
   Data_idx          = np.where(Data_ID_List == idx)[0][0]
   Data_idx_2nd      = np.nonzero(np.isclose(Data_ID_List, idx))[0][0]
   Time_Now          = Time_List[Data_idx_2nd]
   Time_Now_tcross   = Time_List_tcross[Data_idx_2nd]
   if (Data_idx != Data_idx_2nd):
      print("Frame_idx (%i) != np.isclose (%i)"%(Frame_idx, Frame_idx_2nd))
   Frame_idx         = np.where(CM_Time_list == Time_Now)[0][0]
   Frame_idx_2nd     = np.nonzero(np.isclose(CM_Time_list, Time_Now))[0][0]
   Host_CM           = np.array([CM_x_list[Frame_idx_2nd], CM_y_list[Frame_idx_2nd], CM_z_list[Frame_idx_2nd]])
   Host_CM_Rel       = Host_CM - np.array([BOX_SIZE,BOX_SIZE, BOX_SIZE])/2.0
   if (Frame_idx != Frame_idx_2nd):
      print("Frame_idx (%i) != np.isclose (%i)"%(Frame_idx, Frame_idx_2nd))
   NPar_fBound_Now   = float(NPar_Bound[Frame_idx_2nd])/NPar_Bound_Max

   # output z_proj particle density plot
   p = yt.ParticlePlot(ds, ('all','particle_position_x'), ('all','particle_position_y'), ('all','particle_mass'), center='c')
   p.set_font({"size": 25})
   p.set_background_color( field, color="black" )
   p.set_unit('particle_mass', 'code_mass')
   p.hide_colorbar()
   p.hide_axes(draw_frame=True)
   p.annotate_timestamp( time_unit='code_time', corner='upper_left', text_args={'color':'white'}, time_format='[Face-on] $t$ = %.2f crossing time'%Time_Now_tcross)
   p.annotate_text([-(1.174)*(BOX_SIZE/4.0),(1.0)*(BOX_SIZE/4.0)], "$f_\mathrm{bound}$ = %.4f"%NPar_fBound_Now, coord_system="plot", text_args={"color": "white"})
   p.annotate_marker((Host_CM_Rel[0], Host_CM_Rel[1]), coord_system="plot", marker='x', plot_args={"color": "black", "s": 50, "alpha": 0.85})
   p.annotate_sphere(Host_CM, radius=(1, 'code_length'), circle_args={"linewidth": 2.5, "color": "yellow", "alpha": 0.7})#, "linestyle": "--"
   p.annotate_sphere(Host_CM, radius=(10, 'code_length'), circle_args={"linewidth": 3.0, "color": "cyan", "alpha": 0.4})
   p.annotate_text([Host_CM_Rel[0]+1.07, Host_CM_Rel[1]+1.07], "$r_\mathrm{s}$", coord_system="plot", text_args={"color": "yellow", "alpha": 0.85})
   p.annotate_text([Host_CM_Rel[0]+7.5, Host_CM_Rel[1]+7.5], "$10r_\mathrm{s}$", coord_system="plot", text_args={"color": "cyan", "alpha": 0.75})
   p.set_zlim(('all','particle_mass'), (mass_lim_lower, "code_mass"), (mass_lim_upper, "code_mass"))
   p.zoom(zoom_factor)
   p.save(mpl_kwargs={"dpi":dpi})

   # [Zoom] grid structure
   s = yt.SlicePlot(ds, "z", ("gas", "mass"),  center='c')
   s.set_font({"size": 25})
   s.set_background_color( ("gas", "mass"), color="black" )
   s.set_unit('mass', 'code_mass')
   s.set_cmap(field=("gas", "mass"), cmap="gist_gray")
   s.hide_colorbar()
   s.hide_axes(draw_frame=True)
   s.zoom(zoom_factor)
   s.set_zlim(("gas", "mass"), (1e10*UNIT_M, "code_mass"), (1e10*UNIT_M, "code_mass"))
   s.annotate_timestamp( time_unit='code_time', corner='upper_left', text_args={'color':'white'}, time_format='Zoom factor = %d'%zoom_factor)
   s.annotate_grids(alpha=1.0,linewidth=1.0,cmap="binary", edgecolors="white")
   s.annotate_sphere(Host_CM, radius=(1, 'code_length'), circle_args={"linewidth": 2.5, "color": "yellow", "alpha": 0.7})#, "linestyle": "--"
   s.annotate_sphere(Host_CM, radius=(10, 'code_length'), circle_args={"linewidth": 3.0, "color": "cyan", "alpha": 0.4})
   s.annotate_text([Host_CM_Rel[0]+1.07, Host_CM_Rel[1]+1.07], "$r_\mathrm{s}$", coord_system="plot", text_args={"color": "yellow", "alpha": 0.85})
   s.annotate_text([Host_CM_Rel[0]+7.5, Host_CM_Rel[1]+7.5], "$10r_\mathrm{s}$", coord_system="plot", text_args={"color": "cyan", "alpha": 0.75})
   s.save(mpl_kwargs={"dpi":dpi})

   if (idx == idx_start):
      # [Full View] grid structure
      s = yt.SlicePlot(ds, "x", ("gas", "mass"),  center='c')
      s.set_font({"size": 25})
      s.set_background_color( ("gas", "mass"), color="black" )
      s.set_unit('mass', 'code_mass')
      s.set_cmap(field=("gas", "mass"), cmap="gist_gray")
      s.hide_colorbar()
      s.hide_axes(draw_frame=True)
      s.zoom(1)
      s.set_zlim(("gas", "mass"), (1e10*UNIT_M, "code_mass"), (1e10*UNIT_M, "code_mass"))
      s.annotate_timestamp( time_unit='code_time', corner='upper_left', text_args={'color':'white'}, time_format='Zoom factor = 1')
      s.annotate_grids(alpha=1.0,linewidth=1.0,cmap="binary", edgecolors="white")
      s.annotate_sphere(Host_CM, radius=(1, 'code_length'), circle_args={"linewidth": 2.5, "color": "yellow", "alpha": 0.7})#, "linestyle": "--"
      s.annotate_sphere(Host_CM, radius=(10, 'code_length'), circle_args={"linewidth": 3.0, "color": "cyan", "alpha": 0.4})
      s.save(mpl_kwargs={"dpi":dpi})
