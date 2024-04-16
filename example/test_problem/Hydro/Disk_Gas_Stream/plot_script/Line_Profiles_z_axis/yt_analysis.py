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
import matplotlib.ticker as ticker


#-------------------------------------------------------------------------------------------------------------------------
# user-specified parameters
Stars_zlim_upper     = 6.0e8
Stars_zlim_lower     = 1.0e4
Gas_zlim_upper       = 2.5e10
Gas_zlim_lower       = 0.89e2
Gas_vel_lower        = -315.0
Gas_vel_upper        = 315.0
Gas_den_lower        = 5.0e-1
Gas_den_upper        = 1.25e9
Gas_temp_lower       = 0.6e2
Gas_temp_upper       = 4.5e6
Gas_pre_lower        = 2.0e-18
Gas_pre_upper        = 3.0e-9


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

   # load simulation data with yt
   ds          = yt.load('../../Data_%06d'%idx)#, units_override=units_override)
   #print(ds.field_list)
   #print(ds.derived_field_list)
   #print(ds.parameters)
   BoxSize     = float(ds.domain_width[0])#.in_units("kpc"))
   BoxSize_Half= BoxSize/2.0
   Time_Myr    = float(ds.current_time.in_units("Myr"))

   '''
   plot = yt.LinePlot(ds, [("gas", "velocity_x"),("gas", "velocity_y"),("gas", "velocity_z")], [0.0, BoxSize_Half, BoxSize_Half], [BoxSize, BoxSize_Half, BoxSize_Half], 1024)
   plot.annotate_legend([("gas", "velocity_x"),("gas", "velocity_y"),("gas", "velocity_z")])
   plot.set_x_unit("kpc")
   plot.set_unit(("gas", "velocity_x"), "km/s")
   plot.annotate_title(("gas", "velocity_x"), "Gas velocity along x-axis (midplane): %d Myr"%round(Time_Myr))
   plot.set_ylabel("Velocity (km/s)")
   #plot.ylim(('gas', 'velocity_x'), (Gas_vel_lower, "km/s"), (Gas_vel_upper, "km/s"))
   plot.save(mpl_kwargs={"dpi":150})

   plot = yt.LinePlot(ds, ("gas", "density"), [0.0, BoxSize_Half, BoxSize_Half], [BoxSize, BoxSize_Half, BoxSize_Half], 1024)
   plot.annotate_legend(("gas", "density"))
   plot.set_x_unit("kpc")
   plot.set_unit(("gas", "density"), "Msun/kpc**3")
   plot.annotate_title(("gas", "density"), "Gas density along x-axis (midplane): %d Myr"%round(Time_Myr))
   #plot.set_ylim(('gas', 'density'), (Gas_den_lower, "Msun/kpc**3"), (Gas_den_upper, "Msun/kpc**3"))
   plot.save(mpl_kwargs={"dpi":150})

   plot = yt.LinePlot(ds, ("gas", "temperature"), [0.0, BoxSize_Half, BoxSize_Half], [BoxSize, BoxSize_Half, BoxSize_Half], 1024)
   plot.annotate_legend(("gas", "temperature"))
   plot.set_x_unit("kpc")
   plot.set_unit(("gas", "temperature"), "Kelvin")
   plot.annotate_title(("gas", "temperature"), "Gas temperature along x-axis (midplane): %d Myr"%round(Time_Myr))
   #plot.set_ylim(('gas', 'temperature'), (Gas_temp_lower, "Kelvin"), (Gas_temp_upper, "Kelvin"))
   plot.save(mpl_kwargs={"dpi":150})

   plot = yt.LinePlot(ds, ("gas", "pressure"), [0.0, BoxSize_Half, BoxSize_Half], [BoxSize, BoxSize_Half, BoxSize_Half], 1024)
   plot.annotate_legend(("gas", "pressure"))
   plot.set_x_unit("kpc")
   plot.set_unit(("gas", "pressure"), "erg/cm**3")
   plot.annotate_title(("gas", "pressure"), "Gas pressure along x-axis (midplane): %d Myr"%round(Time_Myr))
   #plot.set_ylim(('gas', 'temperature'), (Gas_pre_lower, "erg/cm**3"), (Gas_pre_upper, "erg/cm**3"))
   plot.save(mpl_kwargs={"dpi":150})
   '''

   plot = yt.LinePlot(ds, [("gas", "velocity_x"),("gas", "velocity_y"),("gas", "velocity_z")], [BoxSize_Half, BoxSize_Half, 0.0], [BoxSize_Half, BoxSize_Half, BoxSize], 1024)
   plot.annotate_legend([("gas", "velocity_x"),("gas", "velocity_y"),("gas", "velocity_z")])
   plot.set_x_unit("kpc")
   plot.set_unit(("gas", "velocity_x"), "km/s")
   plot.annotate_title(("gas", "velocity_x"), "Gas velocity along z-axis (zenith): %d Myr"%round(Time_Myr))
   plot.set_ylabel("Velocity (km/s)")
   plot.save(mpl_kwargs={"dpi":150})

   plot = yt.LinePlot(ds, ("gas", "density"), [BoxSize_Half, BoxSize_Half, 0.0], [BoxSize_Half, BoxSize_Half, BoxSize], 1024)
   plot.annotate_legend(("gas", "density"))
   plot.set_x_unit("kpc")
   plot.set_unit(("gas", "density"), "Msun/kpc**3")
   plot.annotate_title(("gas", "density"), "Gas density along z-axis (zenith): %d Myr"%round(Time_Myr))
   plot.save(mpl_kwargs={"dpi":150})

   plot = yt.LinePlot(ds, ("gas", "temperature"), [BoxSize_Half, BoxSize_Half, 0.0], [BoxSize_Half, BoxSize_Half, BoxSize], 1024)
   plot.annotate_legend(("gas", "temperature"))
   plot.set_x_unit("kpc")
   plot.set_unit(("gas", "temperature"), "Kelvin")
   plot.annotate_title(("gas", "temperature"), "Gas temperature along z-axis (zenith): %d Myr"%round(Time_Myr))
   plot.save(mpl_kwargs={"dpi":150})

   plot = yt.LinePlot(ds, ("gas", "pressure"), [BoxSize_Half, BoxSize_Half, 0.0], [BoxSize_Half, BoxSize_Half, BoxSize], 1024)
   plot.annotate_legend(("gas", "pressure"))
   plot.set_x_unit("kpc")
   plot.set_unit(("gas", "pressure"), "erg/cm**3")
   plot.annotate_title(("gas", "pressure"), "Gas pressure along z-axis (zenith): %d Myr"%round(Time_Myr))
   plot.save(mpl_kwargs={"dpi":150})

   '''
   # output star z_proj particle density plot
   p           = yt.ParticlePlot(ds, ('star','particle_position_x'), ('star','particle_position_y'), ('star','particle_mass'), center='c')
   p.set_background_color( ('star','particle_mass'), color="black" )
   p.set_unit(('star','particle_mass'), 'Msun')
   p.hide_colorbar()
   p.hide_axes(draw_frame=True)
   p.annotate_timestamp( time_unit='Myr', corner='upper_left', text_args={'color':'white'}, time_format='[Face-on: Star] $t$ = {time:.1f} {units}')
   p.annotate_marker((0,0), coord_system="plot", marker='+', plot_args={"color": "white", "s": 500, "alpha": 0.85})
   p.set_zlim(('star','particle_mass'), (Stars_zlim_lower, "Msun"), (Stars_zlim_upper, "Msun"))
   p.save(mpl_kwargs={"dpi":150})

   # output star x_proj particle density plot
   p           = yt.ParticlePlot(ds, ('star','particle_position_y'), ('star','particle_position_z'), ('star','particle_mass'), center='c')
   p.set_background_color( ('star','particle_mass'), color="black" )
   p.set_unit(('star','particle_mass'), 'Msun')
   p.hide_colorbar()
   p.hide_axes(draw_frame=True)
   p.annotate_timestamp( time_unit='Myr', corner='upper_left', text_args={'color':'white'}, time_format='[Edge-on: Star] $t$ = {time:.1f} {units}')
   p.annotate_marker((0,0), coord_system="plot", marker='+', plot_args={"color": "white", "s": 500, "alpha": 0.85})
   p.set_zlim(('star','particle_mass'), (Stars_zlim_lower, "Msun"), (Stars_zlim_upper, "Msun"))
   p.save(mpl_kwargs={"dpi":150})

   # output gas z_proj particle density plot
   p           = yt.ProjectionPlot(ds, 'z', ("gas", "density"), center='c')
   #p.annotate_grids(alpha=0.1,linewidth=0.03)
   p.set_background_color( ('gas', 'density'), color="black" )
   p.set_unit(('gas', 'density'), 'Msun/kpc**2')
   p.hide_colorbar()
   p.hide_axes(draw_frame=True)
   p.annotate_timestamp( time_unit='Myr', corner='upper_left', text_args={'color':'white'}, time_format='[Face-on: Gas] $t$ = {time:.1f} {units}')
   p.annotate_marker((0,0), coord_system="plot", marker='+', plot_args={"color": "white", "s": 500, "alpha": 0.85})
   p.set_zlim(('gas', 'density'), (Gas_zlim_lower, "Msun/kpc**2"), (Gas_zlim_upper, "Msun/kpc**2"))
   p.save(mpl_kwargs={"dpi":150})


   # output gas x_proj particle density plot
   p           = yt.ProjectionPlot(ds, 'x', ("gas", "density"), center='c')
   #p.annotate_grids(alpha=0.1,linewidth=0.03)
   p.set_background_color( ('gas', 'density'), color="black" )
   p.set_unit(('gas', 'density'), 'Msun/kpc**2')
   p.hide_colorbar()
   p.hide_axes(draw_frame=True)
   p.annotate_timestamp( time_unit='Myr', corner='upper_left', text_args={'color':'white'}, time_format='[Edge-on: Gas] $t$ = {time:.1f} {units}')
   p.annotate_marker((0,0), coord_system="plot", marker='+', plot_args={"color": "white", "s": 500, "alpha": 0.85})
   p.set_zlim(('gas', 'density'), (Gas_zlim_lower, "Msun/kpc**2"), (Gas_zlim_upper, "Msun/kpc**2"))
   p.save(mpl_kwargs={"dpi":150})
   '''
