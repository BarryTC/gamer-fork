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
Gas_Component        = 1       # [0] off; [1] on
zoom_factor          = 2.5
DM_zlim_upper        = 3.5e8
DM_zlim_lower        = 2.0e6
Stars_zlim_upper     = 6.0e8
Stars_zlim_lower     = 1.0e4
Gas_zlim_upper       = 2.5e10
Gas_zlim_lower       = 0.89e2


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

   def DM_filter(pfilter, data):
      return data[(pfilter.filtered_type, 'ParType')] == 2
   def Stars_filter(pfilter, data):
      return data[(pfilter.filtered_type, 'ParType')] == 3

   yt.add_particle_filter("DM", function=DM_filter, filtered_type="all", requires=["ParType"])
   yt.add_particle_filter("Stars", function=Stars_filter, filtered_type="all", requires=["ParType"])

   # load simulation data with yt
   ds          = yt.load('../../Data_%06d'%idx)#, units_override=units_override)
   BoxSize     = float(ds.domain_width[0].in_units("kpc"))
   #print(BoxSize_Half)
   #sys.exit()
   ds.add_particle_filter("DM")
   ds.add_particle_filter("Stars")
   #print(ds.field_list)
   #print(ds.derived_field_list)

   # [DM] output star z_proj particle density plot
   p           = yt.ParticlePlot(ds, ('DM','particle_position_x'), ('DM','particle_position_y'), ('DM','particle_mass'), center='c')
   p.set_background_color( ('DM','particle_mass'), color="black" )
   p.set_unit(('DM','particle_mass'), 'Msun')
   p.hide_colorbar()
   p.hide_axes(draw_frame=True)
   p.annotate_timestamp( time_unit='Myr', corner='upper_left', text_args={'color':'white'}, time_format='[Face-on: DM] $t$ = {time:.1f} {units}')
   p.annotate_marker([3.0*BoxSize, 3.0*BoxSize], coord_system="plot", marker='+', plot_args={"color": "white", "s": 500, "alpha": 0.85})
   p.zoom(zoom_factor)
   p.set_zlim(('DM','particle_mass'), (DM_zlim_lower, "Msun"), (DM_zlim_upper, "Msun"))
   p.save("DM_%06d_Density_z_Proj.png"%idx, mpl_kwargs={"dpi":150})

   # [DM] output star x_proj particle density plot
   p           = yt.ParticlePlot(ds, ('DM','particle_position_y'), ('DM','particle_position_z'), ('DM','particle_mass'), center='c')
   p.set_background_color( ('DM','particle_mass'), color="black" )
   p.set_unit(('DM','particle_mass'), 'Msun')
   p.hide_colorbar()
   p.hide_axes(draw_frame=True)
   p.annotate_timestamp( time_unit='Myr', corner='upper_left', text_args={'color':'white'}, time_format='[Edge-on: DM] $t$ = {time:.1f} {units}')
   p.annotate_marker([3.0*BoxSize, 3.0*BoxSize], coord_system="plot", marker='+', plot_args={"color": "white", "s": 500, "alpha": 0.85})
   p.zoom(zoom_factor)
   p.set_zlim(('DM','particle_mass'), (DM_zlim_lower, "Msun"), (DM_zlim_upper, "Msun"))
   p.save("DM_%06d_Density_x_Proj.png"%idx, mpl_kwargs={"dpi":150})

   # [Stars] output star z_proj particle density plot
   p           = yt.ParticlePlot(ds, ('Stars','particle_position_x'), ('Stars','particle_position_y'), ('Stars','particle_mass'), center='c')
   p.set_background_color( ('Stars','particle_mass'), color="black" )
   p.set_unit(('Stars','particle_mass'), 'Msun')
   p.hide_colorbar()
   p.hide_axes(draw_frame=True)
   p.annotate_timestamp( time_unit='Myr', corner='upper_left', text_args={'color':'white'}, time_format='[Face-on: Star] $t$ = {time:.1f} {units}')
   p.annotate_marker([3.0*BoxSize, 3.0*BoxSize], coord_system="plot", marker='+', plot_args={"color": "white", "s": 500, "alpha": 0.85})
   p.zoom(zoom_factor)
   p.set_zlim(('Stars','particle_mass'), (Stars_zlim_lower, "Msun"), (Stars_zlim_upper, "Msun"))
   p.save("Stars_%06d_Density_z_Proj.png"%idx, mpl_kwargs={"dpi":150})

   # [Stars] output star x_proj particle density plot
   p           = yt.ParticlePlot(ds, ('Stars','particle_position_y'), ('Stars','particle_position_z'), ('Stars','particle_mass'), center='c')
   p.set_background_color( ('Stars','particle_mass'), color="black" )
   p.set_unit(('Stars','particle_mass'), 'Msun')
   p.hide_colorbar()
   p.hide_axes(draw_frame=True)
   p.annotate_timestamp( time_unit='Myr', corner='upper_left', text_args={'color':'white'}, time_format='[Edge-on: Star] $t$ = {time:.1f} {units}')
   p.annotate_marker([3.0*BoxSize, 3.0*BoxSize], coord_system="plot", marker='+', plot_args={"color": "white", "s": 500, "alpha": 0.85})
   p.zoom(zoom_factor)
   p.set_zlim(('Stars','particle_mass'), (Stars_zlim_lower, "Msun"), (Stars_zlim_upper, "Msun"))
   p.save("Stars_%06d_Density_x_Proj.png"%idx, mpl_kwargs={"dpi":150})

   # [Gas]
   if (Gas_Component == 1):
      # output gas z_proj particle density plot
      p           = yt.ProjectionPlot(ds, 'z', ("gas", "density"), center='c')
      #p.annotate_grids(alpha=0.1,linewidth=0.03)
      p.set_background_color( ('gas', 'density'), color="black" )
      p.set_unit(('gas', 'density'), 'Msun/kpc**2')
      p.hide_colorbar()
      p.hide_axes(draw_frame=True)
      p.annotate_timestamp( time_unit='Myr', corner='upper_left', text_args={'color':'white'}, time_format='[Face-on: Gas] $t$ = {time:.1f} {units}')
      p.annotate_marker([3.0*BoxSize, 3.0*BoxSize], coord_system="plot", marker='+', plot_args={"color": "white", "s": 500, "alpha": 0.85})
      p.zoom(zoom_factor)
      p.set_zlim(('gas', 'density'), (Gas_zlim_lower, "Msun/kpc**2"), (Gas_zlim_upper, "Msun/kpc**2"))
      p.save("Gas_%06d_Density_z_Proj.png"%idx, mpl_kwargs={"dpi":150})

      # output gas x_proj particle density plot
      p           = yt.ProjectionPlot(ds, 'x', ("gas", "density"), center='c')
      #p.annotate_grids(alpha=0.1,linewidth=0.03)
      p.set_background_color( ('gas', 'density'), color="black" )
      p.set_unit(('gas', 'density'), 'Msun/kpc**2')
      p.hide_colorbar()
      p.hide_axes(draw_frame=True)
      p.annotate_timestamp( time_unit='Myr', corner='upper_left', text_args={'color':'white'}, time_format='[Edge-on: Gas] $t$ = {time:.1f} {units}')
      p.annotate_marker([3.0*BoxSize, 3.0*BoxSize], coord_system="plot", marker='+', plot_args={"color": "white", "s": 500, "alpha": 0.85})
      p.zoom(zoom_factor)
      p.set_zlim(('gas', 'density'), (Gas_zlim_lower, "Msun/kpc**2"), (Gas_zlim_upper, "Msun/kpc**2"))
      p.save("Gas_%06d_Density_x_Proj.png"%idx, mpl_kwargs={"dpi":150})
