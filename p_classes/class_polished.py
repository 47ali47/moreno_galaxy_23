import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import copy
from matplotlib import rc
import sys
plt.rcParams['font.family'] = 'serif'


class component_graphs:
 
    def __init__(self):
        
        #container for searching
        self.h = 0.7
        self.particle_data = {}
        self._read_data()
        
################### PRIVATE FUNCTIONS #######################

    def _read_data(self):
    
        read_path = '/Users/anboli/Documents/marckie/object_3094.hdf5'
        
        f = h5py.File(read_path,'r')
    
        #Retrieving specific data and putting it into a container
        
        self.particle_data['dm_x'] = np.asarray(f['dm_x'])
        self.particle_data['dm_y'] = np.asarray(f['dm_y'])
        self.particle_data['dm_z'] = np.asarray(f['dm_z'])
        self.particle_data['dm_mass'] = np.asarray(f['dm_mass'])*1e10/self.h
       
        self.particle_data['gas_x'] = np.asarray(f['gas_x'])
        self.particle_data['gas_y'] = np.asarray(f['gas_y'])
        self.particle_data['gas_z'] = np.asarray(f['gas_z'])
        self.particle_data['gas_mass'] = np.asarray(f['gas_mass'])*1e10/self.h
        
        self.particle_data['stars_x'] = np.asarray(f['stars_x'])
        self.particle_data['stars_y'] = np.asarray(f['stars_y'])
        self.particle_data['stars_z'] = np.asarray(f['stars_z'])
        self.particle_data['stars_mass'] = np.asarray(f['stars_mass'])*1e10/self.h
        self.particle_data['object_Rvir'] = np.asarray(f['object_Rvir'])
        self.particle_data['time'] = np.asarray(f['time'])
        
        f.close()
  
#######################Public_Functions#########################

#this will set the input of plots to whatever variables are put in to var_1 and 2 respectively from the client
    
    def data_arrange(self, component, var_1_name, var_2_name, color, symbol, size, weight):
        
        #setting up the input variables
        var_1 = self.particle_data[component+'_'+var_1_name]
        var_2 = self.particle_data[component+'_'+var_2_name]
        weights = self.particle_data[component+'_'+weight]
    
        
        #code for the graph
        norm = colors.LogNorm(vmin=1e4,vmax=1e9)
        bins = 1024
        cmap = 'nipy_spectral'
        cmap = copy.copy(mpl.cm.get_cmap(cmap))
        hist, x_edges, y_edges = np.histogram2d(var_1, var_2, bins = bins, weights = weights)
        #code to plot
        FOV = 18
        length = 2.*FOV/bins
        area = length**2
        hist = hist/area
        Units = '$M_\u2609$'+' '+ r'$kpc^{-2}$'
        
        plt.imshow(hist, origin = 'lower', extent = [x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]], aspect = 'auto', cmap = cmap, norm = norm)
        hist[hist==0]=1
        cmap.set_bad('black')
        #plotting the circles
        
        
        #inner circle
        theta = np.linspace(0,2 * np.pi, 100)
        radius = self.radius_inner
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)
        plt.plot(x,y)
        
        
        
        #outer circle
        theta_1 = np.linspace(0,2 * np.pi, 100)
        radius = self.radius_outer
        x = radius * np.cos(theta_1)
        y = radius * np.sin(theta_1)
        plt.plot(x,y)
        
       
        #why does this not work here, but works when i put it at the top?
        #plt.rcParams['font.family'] = 'serif'
        rc('text', usetex = True)
        #^^why this only apply to the colorbar
        
        #renaming components such that it appears capitalized in the graph
        if component == 'dm':
            title_name = 'Dark Matter'
        if component == 'stars':
            title_name = 'Stars'
        if component == 'gas':
            title_name = 'Gas'
            
        
        plt.title(f'2D Histogram of {title_name}')
        plt.xlabel(var_1_name + ' '+ '(kpc)')
        plt.ylabel(var_2_name + ' '+ '(kpc)')
        cbar = plt.colorbar()
        #cbar.set_label(r'$\Sigma_{\star}$['+ str(Units) + ']')
        cbar.set_label(r'$\Sigma_{\star}$[$M_{\odot}$ $kpc^{-2}$]')
        plt.savefig(str(component)+"_particle_data.pdf", format="pdf")
        print("open "+str(component)+"_particle_data.pdf")

        