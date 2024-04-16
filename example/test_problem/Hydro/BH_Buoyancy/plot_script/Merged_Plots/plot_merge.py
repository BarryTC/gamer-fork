import argparse
import sys
import yt
import linecache
import numpy as np
import math
import os
import os.path
import matplotlib.image
import matplotlib.pyplot as plt
from numpy import loadtxt


#-------------------------------------------------------------------------------------------------------------------------
# user-specified parameters
folder_path_1   = "../yt_Projection_Plots/"


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
for idx in range(idx_start,idx_end+1,didx):

   print("Data_%06d"%idx)
   # merge subplots
   pz_dens_1 = matplotlib.image.imread(folder_path_1+"Proj_Combine_%06d.png"%idx)
   pz_dens_2 = matplotlib.image.imread("./Data_%06d_Density_Anisotropy_Profiles.png"%idx)
   xdim_1,ydim_1,zdim_1 = pz_dens_1.shape
   xdim_2,ydim_2,zdim_2 = pz_dens_2.shape
   pz_dens_1_zero = np.zeros((max(xdim_1,xdim_2),max(ydim_1,ydim_2),max(zdim_1,zdim_2)))
   pz_dens_2_zero = np.zeros((max(xdim_1,xdim_2),max(ydim_1,ydim_2),max(zdim_1,zdim_2)))
   pz_dens_1_zero[0:xdim_1, 0:ydim_1, 0:zdim_1] = pz_dens_1
   pz_dens_2_zero[0:xdim_2, 0:ydim_2, 0:zdim_2] = pz_dens_2
   merged_image = np.concatenate((pz_dens_1_zero, pz_dens_2_zero), axis=1)
   merged_image = merged_image[0:max(xdim_1,xdim_2),0:ydim_1+ydim_2,:]
   matplotlib.image.imsave("Proj_Combine_%06d.png"%idx, merged_image)
