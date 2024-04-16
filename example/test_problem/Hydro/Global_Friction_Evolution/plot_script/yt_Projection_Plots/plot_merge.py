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
plot1_prefix    = "Data_"
plot1_suffix    = "_Particle_z_particle_mass.png"
plot2_prefix    = "Data_"
plot2_suffix    = "_Slice_z_mass.png"


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
# merge subplots
for idx in range(idx_start,idx_end+1,didx):

    print("Processing snapshot %d/%d"%(idx, idx_end))
    pz_dens_1 = matplotlib.image.imread(plot1_prefix+"%06d"%idx+plot1_suffix)
    pz_dens_2 = matplotlib.image.imread(plot2_prefix+"%06d"%idx+plot2_suffix)
    xdim_1,ydim_1,zdim_1 = pz_dens_1.shape
    xdim_2,ydim_2,zdim_2 = pz_dens_2.shape
    pz_dens_1_zero = np.zeros((max(xdim_1,xdim_2),max(ydim_1,ydim_2),max(zdim_1,zdim_2)))
    pz_dens_2_zero = np.zeros((max(xdim_1,xdim_2),max(ydim_1,ydim_2),max(zdim_1,zdim_2)))
    pz_dens_1_zero[1:xdim_1-1, 1:ydim_1-1, 0:zdim_1] = pz_dens_1[1:xdim_1-1, 1:ydim_1-1, 0:zdim_1]
    pz_dens_2_zero[1:xdim_2-1, 1:ydim_2-1, 0:zdim_2] = pz_dens_2[1:xdim_2-1, 1:ydim_2-1, 0:zdim_2]
    merged_image = np.concatenate((pz_dens_1_zero, pz_dens_2_zero), axis=1)
    matplotlib.image.imsave("Proj_Combine_%06d.png"%idx, merged_image)
