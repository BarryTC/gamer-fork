#!/usr/bin/env python3.9

import argparse
import sys
import numpy as np
import scipy.integrate as integrate
from datetime import datetime
from numpy import loadtxt
import struct


# user-specified parameters
NPar            = int(1e7)
Template        = int(2)                # [1] includes up to ParType; [2] includes ParTag, ParEnergy, ParAngMom
Entry_dtype     = 'f'                   # e.g. 'i', 'f', 'float64'


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

    Input_filename  = "../Particle_%06d.cbin"%idx
    Output_filename = "Particle_%06d.txt"%idx

    # arrays to be outputted
    pmass_list = []
    posx_list  = []
    posy_list  = []
    posz_list  = []
    velx_list  = []
    vely_list  = []
    velz_list  = []
    ptype_list = []
    if (Template == 2):
        ptag_list = []
        peng_list = []
        pam_list  = []


    inh = open(Input_filename, 'rb')
    counter_idx = 0
    while (counter_idx < NPar):
        rec = inh.read(struct.calcsize(Entry_dtype))
        (entry,) = struct.unpack(Entry_dtype, rec)
        pmass_list  += [entry]
        counter_idx += 1
    counter_idx = 0
    while (counter_idx < NPar):
        rec = inh.read(struct.calcsize(Entry_dtype))
        (entry,) = struct.unpack(Entry_dtype, rec)
        posx_list   += [entry]
        counter_idx += 1
    counter_idx = 0
    while (counter_idx < NPar):
        rec = inh.read(struct.calcsize(Entry_dtype))
        (entry,) = struct.unpack(Entry_dtype, rec)
        posy_list   += [entry]
        counter_idx += 1
    counter_idx = 0
    while (counter_idx < NPar):
        rec = inh.read(struct.calcsize(Entry_dtype))
        (entry,) = struct.unpack(Entry_dtype, rec)
        posz_list   += [entry]
        counter_idx += 1
    counter_idx = 0
    while (counter_idx < NPar):
        rec = inh.read(struct.calcsize(Entry_dtype))
        (entry,) = struct.unpack(Entry_dtype, rec)
        velx_list   += [entry]
        counter_idx += 1
    counter_idx = 0
    while (counter_idx < NPar):
        rec = inh.read(struct.calcsize(Entry_dtype))
        (entry,) = struct.unpack(Entry_dtype, rec)
        vely_list   += [entry]
        counter_idx += 1
    counter_idx = 0
    while (counter_idx < NPar):
        rec = inh.read(struct.calcsize(Entry_dtype))
        (entry,) = struct.unpack(Entry_dtype, rec)
        velz_list   += [entry]
        counter_idx += 1
    counter_idx = 0
    while (counter_idx < NPar):
        rec = inh.read(struct.calcsize(Entry_dtype))
        (entry,) = struct.unpack(Entry_dtype, rec)
        ptype_list  += [entry]
        counter_idx += 1
    if (Template == 2):
        counter_idx = 0
        while (counter_idx < NPar):
            rec = inh.read(struct.calcsize(Entry_dtype))
            (entry,) = struct.unpack(Entry_dtype, rec)
            ptag_list   += [entry]
            counter_idx += 1
        counter_idx = 0
        while (counter_idx < NPar):
            rec = inh.read(struct.calcsize(Entry_dtype))
            (entry,) = struct.unpack(Entry_dtype, rec)
            peng_list   += [entry]
            counter_idx += 1
        counter_idx = 0
        while (counter_idx < NPar):
            rec = inh.read(struct.calcsize(Entry_dtype))
            (entry,) = struct.unpack(Entry_dtype, rec)
            pam_list    += [entry]
            counter_idx += 1
    inh.close()


    if (Template == 1):
        with open(Output_filename, 'w') as master_output:
            # output header
            master_output.write("pmass [Msun]             posx [pc]                posy [pc]                posz [pc]                velx [pc/Myr]            vely [pc/Myr]            velz [pc/Myr]            ptype                    \n")
            # output simulation data
            for output_index in range(len(posx_list)):
                    master_output.write(format("%.14e"%pmass_list[output_index],"<25"))
                    master_output.write(format("%.14e"%posx_list[output_index],"<25"))
                    master_output.write(format("%.14e"%posy_list[output_index],"<25"))
                    master_output.write(format("%.14e"%posz_list[output_index],"<25"))
                    master_output.write(format("%.14e"%velx_list[output_index],"<25"))
                    master_output.write(format("%.14e"%vely_list[output_index],"<25"))
                    master_output.write(format("%.14e"%velz_list[output_index],"<25"))
                    master_output.write(format("%.14e"%ptype_list[output_index],"<25"))
                    master_output.write("\n")
    elif (Template == 2):
        with open(Output_filename, 'w') as master_output:
            # output header
            master_output.write("pmass [Msun]             posx [pc]                posy [pc]                posz [pc]                velx [pc/Myr]            vely [pc/Myr]            velz [pc/Myr]            ptype                    ptag                     peng                     pangmom                  \n")
            # output simulation data
            for output_index in range(len(posx_list)):
                    master_output.write(format("%.14e"%pmass_list[output_index],"<25"))
                    master_output.write(format("%.14e"%posx_list[output_index],"<25"))
                    master_output.write(format("%.14e"%posy_list[output_index],"<25"))
                    master_output.write(format("%.14e"%posz_list[output_index],"<25"))
                    master_output.write(format("%.14e"%velx_list[output_index],"<25"))
                    master_output.write(format("%.14e"%vely_list[output_index],"<25"))
                    master_output.write(format("%.14e"%velz_list[output_index],"<25"))
                    master_output.write(format("%.14e"%ptype_list[output_index],"<25"))
                    master_output.write(format("%.14e"%ptag_list[output_index],"<25"))
                    master_output.write(format("%.14e"%peng_list[output_index],"<25"))
                    master_output.write(format("%.14e"%pam_list[output_index],"<25"))
                    master_output.write("\n")