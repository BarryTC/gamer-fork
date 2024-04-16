import sys
import yt
import linecache
import numpy as np
import math
import os
from datetime import datetime
from numpy import loadtxt
print(str(datetime.now()))


#-------------------------------------------------------------------------------------------------------------------------
# user-specified parameters
input_file_path         = "./Omry_Disk_IC/"
input_file_name_list    = ["dm", "stars", "rotcurve"]


#-----------------------------------------------------------------------------------------------------------------
# initial units in 1 = 1e9 Msol/kpc^3 = G = kpc
U_G      = 1.0
U_kpc    = 1.0
U_Msun   = 1.0e-9*(U_kpc**3.0)
# unit conversions
U_pc     = 1.0e-3*U_kpc
U_Mpc    = 1.0e6*U_pc
U_cm     = U_pc /(3.08567758149*1e18)
U_meter  = U_pc /(3.08567758149*1e16)
U_km     = 1.0e3*U_meter
U_gram   = U_Msun/(1.9885*1e33)
U_kg     = U_Msun/(1.9885*1e30)
U_sec    = (6.6738*1e-8*U_cm**3/(U_G*U_gram))**0.5 # gravitational constant
U_yr     = 31556925.2*U_sec                       # tropical year
U_Myr    = 1.0e6*U_yr
U_Gyr    = 1.0e9*U_yr
U_c      = 299792458*U_meter/U_sec               # speed of light
U_eV     = 1.782661907*1e-36*U_kg*U_c**2         # electronvolt
U_MeV    = 1.0e6*U_eV
U_hbar   = 6.582119514*1e-22*U_MeV*U_sec

'''
# constants and conversions in astrophysical units {Msun, pc, Myr} = {1, 1, 1};
U_V_conv = (U_pc/U_Myr)/(U_km/U_sec)
'''


#-------------------------------------------------------------------------------------------------------------------------
# return the numerical value (in str type) of "Parameter" in "FineName"
def Input_Parameter_Readline(FineName, Parameter, extraction_index):
   FileOpen = open('%s'%FineName, 'r').readlines()
   for Line_Now in FileOpen:
      if Parameter in Line_Now:
         return Line_Now.split()[extraction_index]

# compute the Euclidean distance between two spatial points
def distance_cal(reference_origin, point_of_interest):
   origin_x = float(reference_origin[0])
   origin_y = float(reference_origin[1])
   origin_z = float(reference_origin[2])
   pos_x    = float(point_of_interest[0])
   pos_y    = float(point_of_interest[1])
   pos_z    = float(point_of_interest[2])
   return math.sqrt((origin_x-pos_x)**2+(origin_y-pos_y)**2+(origin_z-pos_z)**2)


#-------------------------------------------------------------------------------------------------------------------------
for input_filename_now in input_file_name_list:
   if (input_filename_now == "dm"):
      # load CM position and velocity
      input_data         = loadtxt(input_file_path+input_filename_now, dtype=float)
      Pos_x_list         = input_data[:,0].astype(float)
      Pos_y_list         = input_data[:,1].astype(float)
      Pos_z_list         = input_data[:,2].astype(float)
      Vel_x_list         = input_data[:,3].astype(float)
      Vel_y_list         = input_data[:,4].astype(float)
      Vel_z_list         = input_data[:,5].astype(float)
      Mass_list          = input_data[:,6].astype(float)

      # output processed data table
      '''
      processed_r_kpc = Par_List_Bin_R*1e-3
      processed_beta  = list(zip(*Par_Attr_List_Proc))[0]
      processed_v3D   = list(zip(*Par_Attr_List_Proc))[2]
      processed_rho   = Par_List_Bin_rho*1e9
      processed_par_n = list(zip(*Par_Attr_List_Proc))[1]
      '''
      ListLen            = len(Mass_list)
      with open("DM_IC_List" , 'w') as writer:
         # output header
         writer.write("Radius [kpc]             Pos_x [pc]               Pos_y [pc]               Pos_z [pc]               Vel_x [km/sec]           Vel_y [km/sec]           Vel_z [km/sec]           \n")
         # output simulation data
         for output_index in range(ListLen):
            writer.write(format("%.14e"%Mass_list[output_index],"<25"))
            writer.write(format("%.14e"%Pos_x_list[output_index],"<25"))
            writer.write(format("%.14e"%Pos_y_list[output_index],"<25"))
            writer.write(format("%.14e"%Pos_z_list[output_index],"<25"))
            writer.write(format("%.14e"%Vel_x_list[output_index],"<25"))
            writer.write(format("%.14e"%Vel_y_list[output_index],"<25"))
            writer.write(format("%.14e"%Vel_z_list[output_index],"<25"))
            writer.write("\n")
      print("DM_IC_List complete: N = %d"%ListLen)
      print(str(datetime.now()))
   elif (input_filename_now == "stars"):
      # load CM position and velocity
      input_data         = loadtxt(input_file_path+input_filename_now, dtype=float)
      Pos_x_list         = input_data[:,0].astype(float)
      Pos_y_list         = input_data[:,1].astype(float)
      Pos_z_list         = input_data[:,2].astype(float)
      Vel_x_list         = input_data[:,3].astype(float)
      Vel_y_list         = input_data[:,4].astype(float)
      Vel_z_list         = input_data[:,5].astype(float)
      Mass_list          = input_data[:,6].astype(float)

      # output processed data table
      '''
      processed_r_kpc = Par_List_Bin_R*1e-3
      processed_beta  = list(zip(*Par_Attr_List_Proc))[0]
      processed_v3D   = list(zip(*Par_Attr_List_Proc))[2]
      processed_rho   = Par_List_Bin_rho*1e9
      processed_par_n = list(zip(*Par_Attr_List_Proc))[1]
      '''
      ListLen            = len(Mass_list)
      with open("Stars_IC_List" , 'w') as writer:
         # output header
         writer.write("Radius [kpc]             Pos_x [pc]               Pos_y [pc]               Pos_z [pc]               Vel_x [km/sec]           Vel_y [km/sec]           Vel_z [km/sec]           \n")
         # output simulation data
         for output_index in range(ListLen):
            writer.write(format("%.14e"%Mass_list[output_index],"<25"))
            writer.write(format("%.14e"%Pos_x_list[output_index],"<25"))
            writer.write(format("%.14e"%Pos_y_list[output_index],"<25"))
            writer.write(format("%.14e"%Pos_z_list[output_index],"<25"))
            writer.write(format("%.14e"%Vel_x_list[output_index],"<25"))
            writer.write(format("%.14e"%Vel_y_list[output_index],"<25"))
            writer.write(format("%.14e"%Vel_z_list[output_index],"<25"))
            writer.write("\n")
      print("Stars_IC_List complete: N = %d"%ListLen)
      print(str(datetime.now()))
   elif (input_filename_now == "rotcurve"):
      # load CM position and velocity
      input_data         = loadtxt(input_file_path+input_filename_now, dtype=float)
      Radius_list        = input_data[:,0].astype(float)
      Rot_Vel_list       = input_data[:,1].astype(float)

      # output processed data table
      '''
      processed_r_kpc = Par_List_Bin_R*1e-3
      processed_beta  = list(zip(*Par_Attr_List_Proc))[0]
      processed_v3D   = list(zip(*Par_Attr_List_Proc))[2]
      processed_rho   = Par_List_Bin_rho*1e9
      processed_par_n = list(zip(*Par_Attr_List_Proc))[1]
      '''
      ListLen            = len(Radius_list)
      with open("RotCurve_IC_List" , 'w') as writer:
         # output header
         writer.write("Radius [kpc]             Rot_Vel [km/sec]         \n")
         # output simulation data
         for output_index in range(ListLen):
            writer.write(format("%.14e"%Radius_list[output_index],"<25"))
            writer.write(format("%.14e"%Rot_Vel_list[output_index],"<25"))
            writer.write("\n")
      print("RotCurve_IC_List complete: N = %d"%ListLen)
      print(str(datetime.now()))
