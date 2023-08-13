
import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import copy
from matplotlib import rc
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
        

    def total_mass(self):
        
        s_mass  = self.particle_data['stars_mass']
        dm_mass = self.particle_data['dm_mass']
        g_mass  = self.particle_data['gas_mass']
        
        s_total_mass  = sum(s_mass)
        dm_total_mass = sum(dm_mass)
        g_total_mass  = sum(g_mass)
        
        #print('s total mass =', s_total_mass)
        #print('dm total mass =', dm_total_mass)
        #print('g total mass =', g_total_mass)
    
    #def star_mask(self):
        #for

     
    def mask(self, radius_inner, radius_outer, component):

        self.radius_inner = radius_inner
        self.radius_outer = radius_outer
        
        x_input = component + '_x'
        y_input = component + '_y'
        r_input = component + '_r'
        input_dictionary = {}
       
        #Here I am setting a key value pair to the component variable with its corresponding list of values
        
        input_dictionary[x_input] = self.particle_data[component+'_x']
        input_dictionary[y_input] = self.particle_data[component+'_y']
        input_dictionary[r_input] = np.sqrt(input_dictionary[x_input]**2+input_dictionary[y_input]**2)
        
        #the values of the radius of each point is stored in this dictionary
        print('input dictionary r-input')
        print((input_dictionary[r_input]))
        
        #testing
        print('length of stars-mass')
        print(len(self.particle_data['stars_mass']))
        print('len of input dictionary r-input')
        print(len(input_dictionary[r_input]))
        
        #defining mask to be the r values between r_inner and r_outer (rings)
        
        self.mask = (self.radius_inner < input_dictionary[r_input]) & (self.radius_outer > input_dictionary[r_input])
        #print(self.mask)
        
        
        #print('hello')
        #print((self.particle_data['stars_mass']))
    
    def particles_ring(self, component):
        
        self.component = component
        m_variable = 'm_' + component
        m_dictionary = {}
        #assigning key(m_component) to corresponding mask of particles within the radius determined above
        m_dictionary[m_variable] = np.sum(self.mask)
        #print(f'sum of {component} between radius {self.radius_inner} and {self.radius_outer} = '+ str(m_dictionary[m_variable]))
        
    def star_mass_ring(self, component):
        #sums up the mass of stars between rings
        print(sum(self.particle_data['stars_mass'][self.mask]))
        
        
    def star_mass_radius_plot(self):
        x_input = 'stars_x'
        y_input = 'stars_y'
        r_input = 'stars_r'
        
        input_dictionary = {}
       
        #Here I am setting a key value pair to the component variable with its corresponding list of values
        
        input_dictionary[x_input] = self.particle_data['stars_x']
        input_dictionary[y_input] = self.particle_data['stars_y']
        input_dictionary[r_input] = np.sqrt(input_dictionary[x_input]**2+input_dictionary[y_input]**2)
        
        
        
        #dictionary to store radius key to mass value
        radius_mass = {}
        print('joey')
        for i in range(36043):
            radius_mass[input_dictionary[r_input][i]] = self.particle_data['stars_mass'][i]
        
        keys = list(radius_mass.keys())
        values = tuple(list(radius_mass.values()))
        
        x_values = [pair[0] for pair in values]
        y_values = [pair[1] for pair in values]
        plt.figure(figsize=(10, 6))  # Set the size of the figure

        # Create the scatter plot
        plt.scatter(x_values, y_values, color='blue', label='Data Points')

        # Customize the plot
        plt.title('Scatter Plot of Radius vs. Mass')
        plt.xlabel('Radius')
        plt.ylabel('Mass')
        plt.legend()

        # Display the plot
        plt.show()
         
         
         
        #print(input_dictionary[r_input][36042])
        #print(self.particle_data['stars_mass'][36042])
        
        
        
            
        
        

        
       
        #print('hello1')
        #print(np.sum(self.mask))
        #print('hello2')
        #print(self.particle_data['stars_mass'][self.mask])
        #print(sum(self.particle_data['stars_mass'][self.mask]))
        #print('length is',len(self.particle_data['stars_mass'][self.mask]) )
        #print(len(self.particle_data['stars_mass'][self.mask]))
        #checking to see the length of the data set checks out :)
        #print(len(self.particle_data['stars_mass']))
        #print(len(self.particle_data['stars_x']))
        #print(len(self.mask))
        
        
        
        
        
   

        
        
        
        
        
        
       
"""
    def mask(self, radius_inner, radius_outer, component):

        self.radius_inner = radius_inner
        self.radius_outer = radius_outer
        
        x_input = component + '_x'
        y_input = component + '_y'
        r_input = component + '_r'
        x_value = self.particle_data[component+'_x']
        y_value = self.particle_data[component+'_y']
        r_input_value = np.sqrt(x_value**2+y_value**2)
        
       
        #Here I am setting a key value pair to the component variable with its corresponding list of values
        
        input_dictionary = {}
        input_dictionary[x_input] = x_value
        input_dictionary[y_input] = y_value
        input_dictionary[r_input] = r_input_value
        
        #the values of the radius of each point is stored in this dictionary
        print(input_dictionary[r_input])
        
        #defining mask to be the r values between r_inner and r_outer (rings)
        mask = (self.radius_inner < input_dictionary[r_input]) & (self.radius_outer > input_dictionary[r_input])
        print(mask)
"""
        
        
        
        
        
        
        
        
        
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
