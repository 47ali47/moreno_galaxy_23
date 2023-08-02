
import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

##################### CLASS DEFINITION ######################

class Human:
 
    def __init__(self):
        
        #container for searching
        self.cosmic_data = {}
        self._read_data()
        #self._plot_data()
        
################### PRIVATE FUNCTIONS #######################

    def _read_data(self):
    
        read_path = '/Users/anboli/Documents/marckie/global_sample_data_snapshot_1200.hdf5'
        
        f = h5py.File(read_path,'r')
    
        #Retrieving specific data (x, y, radius, rvir) and putting it into a container
        self.cosmic_data['Rvir'] = np.asarray(f['Rvir'])
        self.cosmic_data['x'] = np.asarray(f['Xc_ahf_cat'])
        self.cosmic_data['y'] = np.asarray(f['Yc_ahf_cat'])
        self.cosmic_data['r'] = np.asarray(f['r80_stars'])
        self.cosmic_data['z'] = np.asarray(f['Zc_ahf_cat'])
        self.cosmic_data['galaxyID'] = np.asarray(f['galaxyID'])
        self.cosmic_data['groupID'] = np.asarray(f['groupID'])
        
        #Here, i is the index, rvir, and the x and y also get printed at each specific index i
        
        #temp turned off
        #for i in np.arange(self.cosmic_data['Rvir'].shape[0]):
            #print(i,self.cosmic_data['Rvir'][i],self.cosmic_data['x'][i],self.cosmic_data['y'][i],self.cosmic_data['z'][i] )
        
        print(self.cosmic_data['galaxyID'].shape[0])
        print(self.cosmic_data['galaxyID'].shape)
        
        
        for i in np.arange(self.cosmic_data['galaxyID'].shape[0]):
            galaxyID = self.cosmic_data['galaxyID'][i]
            groupID  = self.cosmic_data['groupID'][i]
            if groupID == -1.0:
                print(i,
                    galaxyID,
                    groupID,
                    groupID == -1.0)
                  
        
        
        
        #for i, value in enumerate(data):
            #if mask[i]:
                #print(i, value)
    
        #for i in mask_true:
            #print(i, self.cosmic_data['galaxyID'][i], self.cosmic_data['groupID'][i],mask_true[i])
        
        
            #if mask == True:
                
          #original copy
#        for i in np.arange(self.cosmic_data['galaxyID'].shape[0]):
#            print(i, self.cosmic_data['galaxyID'][i], self.cosmic_data['groupID'][i], mask_central[i])
        
        #for i in np.arange(self.cosmic_data['galaxyID'].shape[0]):
            #print(i, self.cosmic_data['galaxyID'][i], self.cosmic_data['groupID'][i])
            
#                    if mask == True:
#                for i in np.arange(self.cosmic_data['galaxyID'].shape[0]):
#                print(mask[i])

        
        
        #to store memory
        f.close()
        
########################################################################
        
    #def _plot_data(self):
        
        #renaming variables to make it prettier
       
        #x = self.cosmic_data['x']
        #y = self.cosmic_data['y']
        #z = self.cosmic_data['y']
        #Rvir = self.cosmic_data['Rvir']
        #radius = self.cosmic_data['r']
        
        #scale by 0.7 because of hubble constant
        #plt.scatter(x/0.7, y/0.7, s = 0.1*np.log10(Rvir), c = radius)
       
        #Graph setup need to do z axis
        #plt.title('Coordinates')
        #plt.ylabel('y (kpc)')
        #plt.xlabel('x (kpc)')
        #plt.legend()
        #plt.savefig("marckiedata.pdf", format="pdf")
        #print("open marckiedata.pdf")
        
        #you need to save the figure before you show for easier access

#######################Public_Functions#########################

#this will set the input of plots to whatever variables are put in to var_1 and 2 respectively from the client
    
    def data_arrange(self, var_1_name, var_2_name, color, symbol, size):
    
        Rvir = self.cosmic_data['Rvir']
        radius = self.cosmic_data['r']
        
        #if x is entered, it will search for the x variable etc...
        var_1 = self.cosmic_data[var_1_name]
        var_2 = self.cosmic_data[var_2_name]
        
        plt.scatter(var_1/0.7, var_2/0.7, s = size, c = color, marker = symbol)
        plt.title('Coordinates')
        plt.xlabel(var_1_name + '(kpc)')
        plt.ylabel(var_2_name + '(kpc)')
        plt.legend()
        plt.savefig("marckiedata.pdf", format="pdf")
        print("open marckiedata.pdf")
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
#        self.var_1 = var_1_name
#        self.var_2 = var_2_name
#        self.color = color
#        self.symbol = symbol
#        self.size = size
#        Rvir = self.cosmic_data['Rvir']
#        #radius = self.cosmic_data['r']
#
#
#        print(('var1=')+(var_1_name)+('   ')+('var2=')+(var_2_name)+('   ')+('color=')+(color)+(' ' )+('symbol=')+(symbol)+(' ')+('size=')+(str(size)))
#
#        plt.scatter(float(self.var_1)/0.7,float(self.var_2)/0.7, marker = symbol, s = size, c = color)
#        plt.title('Coordinates')
#        plt.xlabel(var_1_name + '(kpc)')
#        plt.ylabel(var_2_name + '(kpc)')
#        plt.legend()
#        plt.savefig("marckiedata.pdf", format="pdf")
#        print("open marckiedata.pdf")
#
#        #x_variables
#
#        #var_1 = self.cosmic_data[var_1_name]
#        #var_2 = self.cosmic_data[var_2_name]
#"""
#        if var_1 == 'x':
#            var_1_new = self.cosmic_data['x']
#            var_1_name = 'x'
#        if var_1 == 'y':
#            var_1_new = self.cosmic_data['y']
#            var_1_name = 'y'
#        if var_1 == 'z':
#            var_1_new = self.cosmic_data['z']
#            var_1_name = 'z'
#
#        #y_variables
#
#        if var_2 == 'x':
#            var_2_new = self.cosmic_data['x']
#            var_2_name = 'x'
#        if var_2 == 'y':
#            var_2_new = self.cosmic_data['y']
#            var_2_name = 'y'
#        if var_2 == 'z':
#            var_2_new = self.cosmic_data['z']
#            var_2_name = 'z'
#        """
#
#        #just to see if the plugins are correct
#
#
#
##    def plot_maps(self):
##
##        views = ['xy','yz','zx' ]
##        for view in views:
##
##            if view == 'xy':
##                var_1 = 'x'
##                var_2 = 'y'
##                plt.title('Coordinates')
##                plt.xlabel(var_1 + '(kpc)')
##                plt.ylabel(var_2 + '(kpc)')
##                plt.legend()
##                plt.scatter(var_1/0.7, var_2/0.7, marker = symbol, s = size, c = color)
##                plt.savefig("marckiedata_xy.pdf", format="pdf")
##                print("open marckiedata_xy.pdf")
##
##            if view == 'yz':
##                var_1 = 'y'
##                var_2 = 'z'
##                plt.title('Coordinates')
##                plt.xlabel(var_1 + '(kpc)')
##                plt.ylabel(var_2 + '(kpc)')
##                plt.legend()
##                plt.scatter(var_1/0.7, var_2/0.7, marker = symbol, s = size, c = color)
##                plt.savefig("marckiedata_yz.pdf", format="pdf")
##                print("open marckiedata_yz.pdf")
##
##            if view == 'zx':
##                var_1 = 'z'
##                var_2 = 'x'
##                plt.title('Coordinates')
##                plt.xlabel(var_1 + '(kpc)')
##                plt.ylabel(var_2 + '(kpc)')
##                plt.legend()
##                plt.scatter(var_1/0.7, var_2/0.7, marker = symbol, s = size, c = color)
##                plt.savefig("marckiedata_zx.pdf", format="pdf")
##                print("open marckiedata_zx.pdf")
#
#
#
#
#
#
