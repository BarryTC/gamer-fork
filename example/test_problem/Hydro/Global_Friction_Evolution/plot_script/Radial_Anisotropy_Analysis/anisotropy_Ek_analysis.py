import argparse
import sys
import yt
import linecache
import numpy as np
import math
import os
import matplotlib.image
import matplotlib.pyplot as plt
from datetime import datetime
from numpy import loadtxt
from operator import add
print(str(datetime.now()))


#-------------------------------------------------------------------------------------------------------------------------
# user-specified parameters
radial_Nbin      = 50


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

# compute the Euclidean distance between two spatial points
def distance_cal(reference_origin, point_of_interest):
   origin_x = float(reference_origin[0])
   origin_y = float(reference_origin[1])
   origin_z = float(reference_origin[2])
   pos_x    = float(point_of_interest[0])
   pos_y    = float(point_of_interest[1])
   pos_z    = float(point_of_interest[2])
   return math.sqrt((origin_x-pos_x)**2+(origin_y-pos_y)**2+(origin_z-pos_z)**2)

# parse strings of particle attributes
def par_attr_parse(stringline, CM_coords, CM_vel):#, radial_log10_bin_start, radial_log10_bin_width):
   pos_vector = np.array([float(i) for i in stringline.split()[1:4]])-CM_coords
   vel_vector = np.array([float(i) for i in stringline.split()[4:7]])-CM_vel
   # (1) compute distance to the system's CM
   To_CM_distance = distance_cal([0.0,0.0,0.0], pos_vector)
   # (2) determine the radial component of the velocity via vector projection
   vel_rad = pos_vector*np.dot(vel_vector, pos_vector)/np.dot(pos_vector, pos_vector)
   # (3) determine the tangential component of the velocity via vector subtraction
   vel_tan = np.subtract(vel_vector, vel_rad)
   if (np.linalg.norm(vel_tan) > np.linalg.norm(vel_vector)):
      print("Error!!")
   # (4) output [radial distance to CM, \sigma_{radial}^2, \sigma_{tangential}^2]
   output_list = [To_CM_distance, np.linalg.norm(vel_rad)**2, np.linalg.norm(vel_tan)**2]
   return output_list


#-------------------------------------------------------------------------------------------------------------------------
# load problem-specific parameters
Input_TestProb        = "../../PAR_IC__Input__TestProb"
BoxSize_Half          = float(Input_Parameter_Readline(Input_TestProb, "BOX_SIZE    ", 1))/2.0
NFW_concentration     = float(Input_Parameter_Readline(Input_TestProb, "NFW_concentration    ", 1))
Input_TestProb        = "../../Input__TestProb"
ParM                  = float(Input_Parameter_Readline(Input_TestProb, "Particle_Mass    ", 1))

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


#-------------------------------------------------------------------------------------------------------------------------
# begin data processing
for idx in range(idx_start,idx_end+1,didx):
   print(str(datetime.now()))
   print("Start processing Frame %i/%i"%(idx, idx_end))
   # subhalo CM position
   Data_idx           = np.where(Data_ID_List == idx)[0][0]
   Data_idx_2nd       = np.nonzero(np.isclose(Data_ID_List, idx))[0][0]
   Time_Now           = Time_List[Data_idx_2nd]
   if (Data_idx != Data_idx_2nd):
      print("Check1: Frame_idx (%i) != np.isclose (%i)"%(Frame_idx, Frame_idx_2nd))
   Frame_idx          = np.where(CM_Time_list == Time_Now)[0][0]
   Frame_idx_2nd      = np.nonzero(np.isclose(CM_Time_list, Time_Now))[0][0]
   Coord_CoM_Pos_Vec  = np.array([CM_x_list[Frame_idx], CM_y_list[Frame_idx], CM_z_list[Frame_idx]])
   Vel_CoM_Pos_Vec    = np.array([CM_Vx_list[Frame_idx], CM_Vy_list[Frame_idx], CM_Vz_list[Frame_idx]])
   if (Frame_idx != Frame_idx_2nd):
      print("Check2: Frame_idx (%i) != np.isclose (%i)"%(Frame_idx, Frame_idx_2nd))
   NPar_Bound_Now     = int(NPar_Bound_list[Frame_idx_2nd])

   # load the targeted "Particle_xxxxxx.txt" file
   Par_Input_Filename = "../Particle_%06d.txt"%idx
   ParTag             = loadtxt(Par_Input_Filename, skiprows=1, dtype=float)[:,8].astype(float)
   Par_Bound_List     = np.array(range(int(len(ParTag))))[ParTag > 0.0] + 1 + 1
   print("NPar_Bound_Now = %i, len(Par_Bound_List) = %i"%(NPar_Bound_Now, len(Par_Bound_List)))
   NPar_Bound_Now     = len(Par_Bound_List)
   Par_Ind_Attr_List  = []                                                     # [CoM distance, \sigma_{radial}^2, \sigma_{tangential}^2]
   Par_Attr_List      = []                                                     # [particle count, \Sum\sigma_{radial}^2, \Sum\sigma_{tangential}^2]
   Par_List_Bin_R     = np.array([])                                           # record corresponding radius of each bins
   Par_List_Bin_rho   = np.array([])                                           # record averaged density enclosed within each bin

   # select bound particles only
   for par_idx in range(NPar_Bound_Now):
      Par_Ind_Attr_List += [par_attr_parse(linecache.getline(Par_Input_Filename, Par_Bound_List[par_idx]), Coord_CoM_Pos_Vec, Vel_CoM_Pos_Vec)]
   Par_Ind_Attr_List  = np.array(Par_Ind_Attr_List)
   # sort the particle entries based on CoM distance
   Par_Ind_Attr_List  = Par_Ind_Attr_List[np.argsort(Par_Ind_Attr_List[:,0])]

   # start grouping particle entries by the respective CoM distances
   ini_bin_entry_num  = [20, 20, 25, 30, 45, 50]
   initial_bin_len    = len(ini_bin_entry_num)
   cum_elem_count     = 0
   distance_i_shell   = 0
   distance_f_shell   = 0
   # initial five radial bins
   for ini_bin_size in ini_bin_entry_num:
      # [particle count, \Sum\sigma_{radial}^2, \Sum\sigma_{tangential}^2]
      Par_Attr_List += [[0.0, 0.0, 0.0]]
      for idx_per_bin in range(ini_bin_size):
         Par_Attr_List[-1][0] += 1
         Par_Attr_List[-1][1] += Par_Ind_Attr_List[cum_elem_count+idx_per_bin][1]
         Par_Attr_List[-1][2] += Par_Ind_Attr_List[cum_elem_count+idx_per_bin][2]
      # update parameters after the curret iteration
      cum_elem_count  += ini_bin_size
      distance_f_shell = (Par_Ind_Attr_List[cum_elem_count-1][0]+Par_Ind_Attr_List[cum_elem_count][0])/2
      Par_List_Bin_R   = np.append(Par_List_Bin_R, (distance_i_shell+distance_f_shell)/2)
      Par_List_Bin_rho = np.append(Par_List_Bin_rho, ini_bin_size*ParM/(4*math.pi*((distance_f_shell)**3-(distance_i_shell)**3)/3))
      distance_i_shell = distance_f_shell

   # subsequent radial bins
   Radial_Bin_Start  = math.log10(distance_f_shell)
   Radial_Bin_End    = math.log10(BoxSize_Half)
   Radial_Bin_Width  = (Radial_Bin_End-Radial_Bin_Start)/radial_Nbin
   Radial_Bin_log10  = np.linspace(Radial_Bin_Start, Radial_Bin_End, radial_Nbin+1)
   # [particle count, \Sum\sigma_{radial}^2, \Sum\sigma_{tangential}^2]
   Par_Attr_List += [[0.0, 0.0, 0.0]]*(radial_Nbin+1)
   Par_Attr_List  = np.array(Par_Attr_List)
   Par_Attr_len   = np.shape(Par_Attr_List)[0]-1
   for par_idx in range(cum_elem_count, NPar_Bound_Now):
      # register attributes for particles with distances within the targeted range
      bin_idx = int(math.floor((math.log10(Par_Ind_Attr_List[par_idx][0])-Radial_Bin_Start)/Radial_Bin_Width))
      if (bin_idx+initial_bin_len <= Par_Attr_len):
         Par_Attr_List[bin_idx+initial_bin_len][0] += 1
         Par_Attr_List[bin_idx+initial_bin_len][1] += Par_Ind_Attr_List[par_idx][1]
         Par_Attr_List[bin_idx+initial_bin_len][2] += Par_Ind_Attr_List[par_idx][2]
   # update parameters after the curret iteration
   for bin_idx in range(initial_bin_len, Par_Attr_len):
      Par_List_Bin_R   = np.append(Par_List_Bin_R, 10**((Radial_Bin_log10[bin_idx-initial_bin_len]+Radial_Bin_log10[bin_idx-initial_bin_len+1])/2))
      distance_f_shell = 10**Radial_Bin_log10[bin_idx-initial_bin_len+1]
      Par_List_Bin_rho = np.append(Par_List_Bin_rho, Par_Attr_List[bin_idx][0]*ParM/(4*math.pi*((distance_f_shell)**3-(distance_i_shell)**3)/3))
      distance_i_shell = distance_f_shell

   # compute enclosed mass distribution
   M_enc_List          = []
   for bin_idx, R_now in enumerate(Par_List_Bin_R):
      M_enc_List      += [len(Par_Ind_Attr_List[Par_Ind_Attr_List[:,0] < R_now])]
   M_enc_List          = ParM*np.array(M_enc_List)

   # output ensemble-averaged radial anisotropy in each radial bin
   Par_Attr_List_Proc = [[0.0, 0, 0.0] for _ in range(Par_Attr_len)] # [beta, particle count, ave(\sigma_{total}^2)]
   for bin_idx in reversed(range(Par_Attr_len)):
      [par_num_total, sigma_rad_sqr_total, sigma_tan_sqr_total] = Par_Attr_List[bin_idx]
      # remove zero entries
      if (par_num_total < 1):
         Par_List_Bin_R     = np.delete(Par_List_Bin_R, bin_idx)
         Par_List_Bin_rho   = np.delete(Par_List_Bin_rho, bin_idx)
         Par_Attr_List_Proc = np.delete(Par_Attr_List_Proc, obj=bin_idx, axis=0)
         M_enc_List         = np.delete(M_enc_List, bin_idx)
      else:
         Par_Attr_List_Proc[bin_idx][0] = 1- sigma_tan_sqr_total/(2*sigma_rad_sqr_total)
         Par_Attr_List_Proc[bin_idx][1] = par_num_total
         Par_Attr_List_Proc[bin_idx][2] = ((sigma_tan_sqr_total+sigma_rad_sqr_total)/par_num_total)**0.5

   # output processed data table
   processed_r     = Par_List_Bin_R
   processed_beta  = np.array(list(zip(*Par_Attr_List_Proc))[0])
   processed_v3D   = np.array(list(zip(*Par_Attr_List_Proc))[2])
   processed_rho   = Par_List_Bin_rho
   processed_par_n = list(zip(*Par_Attr_List_Proc))[1]
   processed_M_enc = M_enc_List
   with open("Particle_Radial_Anisotropy_%06d"%idx , 'w') as writer:
      # output header
      writer.write("Radius [Code_L]          Anisotropy               3D Velocity [Code_V]     Density [Code_rho]       Particle Counts          Enclosed Mass [Code_M]        \n")
      # output simulation data
      for output_index in range(len(processed_r)):
         writer.write(format("%e"%processed_r[output_index],"<25"))
         writer.write(format("%e"%processed_beta[output_index],"<25"))
         writer.write(format("%e"%processed_v3D[output_index],"<25"))
         writer.write(format("%e"%processed_rho[output_index],"<25"))
         writer.write(format("%d"%processed_par_n[output_index],"<25"))
         writer.write(format("%e"%processed_M_enc[output_index],"<25"))
         writer.write("\n")
   print("Particle_Radial_Anisotropy_%06d complete"%idx)
