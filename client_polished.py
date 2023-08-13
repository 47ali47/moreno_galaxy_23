import math
import p_classes.class_polished as class_polished
import sys
import os
import subprocess
import numpy as np
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import copy
from matplotlib import rc
plt.rcParams['font.family'] = 'serif'
data_obj = class_polished.component_graphs()

stars_x = data_obj.particle_data['stars_x']
stars_y = data_obj.particle_data['stars_y']
stars_mass = data_obj.particle_data['stars_mass']

def radius(x, y):
    return math.sqrt(x**2+y**2)


def calculate_ring_area(inner_radius, outer_radius):
    return math.pi * ( outer_radius**2 - inner_radius**2)


def handle_ring(ring_criteria, all_stars):
    """
    Identify stars within the ring according to ring_criteria. 
    Then return the number of stars in the ring and their total mass and their mass density.
    :param ring_criteria: (inner_radius, outer_radius)
    :param all_stars: [((x, y), mass)...((x, y)), mass)]
    :return: (number_of_stars_in_the_ring, total_star_mass_in_the_ring, ring_mass_density)
    """
    (radius_inner, radius_outer) = ring_criteria
    star_count = 0
    total_star_mass = 0
    ring_area = calculate_ring_area(radius_inner, radius_outer) 
    #identify stars in the ring, star count, and their total mass
    for star in all_stars:
        #unpacking the variables
        ((x,y), mass) = star
        if radius_inner < radius(x,y) <= radius_outer: #inside the ring
            star_count = star_count + 1
            total_star_mass = total_star_mass + mass
    return (star_count, total_star_mass, total_star_mass/ring_area)

#all_stars is a list of the stars coordinate and mass
#[((x, y), mass)...((x, y)), mass)]
all_stars = []

rings = {} #{ring_index:(star_count, total_star_mass_in_ring, ring_density)}

#length is 36043 the number of the stars
for i in range(len(stars_x)):
    coordinate = (stars_x[i],stars_y[i])
    mass = stars_mass[i]
    all_stars.append((coordinate, mass))

#max radius = 17.96
max_radius = data_obj.particle_data['object_Rvir']
slice_thickness = 0.1
number_slices = math.ceil(max_radius/slice_thickness)

#defining the rings
radius_inner = 0
radius_outer = slice_thickness
for slice_index in range(number_slices):
    rings[slice_index] = handle_ring((radius_inner, radius_outer), all_stars)
    radius_inner = radius_outer
    radius_outer = radius_inner + slice_thickness


def plot_star_rings(filename, title, x_label, x_values, y_label, y_values, label,
                    color = 'blue', figsize = (10,6)):
    plt.figure(figsize=figsize)  # Set the size of the figure

    # Create the scatter plot
    plt.scatter(x_values, y_values, color=color, label=label)

    # Customize the plot
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend()

    # Display the plot
    plt.savefig(filename, format="pdf")
    os.system("open {}".format(filename))

#StarMass graph
plot_label = 'Ring Thickness = {}, \n Max Radius = {} KPC'.format(slice_thickness, max_radius)
x_values = rings.keys() #the ring indexes e.g. 0, 1, 2, 3...
y_values = [mass for (count, mass, density) in rings.values()]
plot_star_rings(filename = 'star_mass_plot.pdf',
                title = 'Scatter Plot of Star Radius (in rings) vs. Mass',
                x_label = 'Rings',
                x_values = x_values,
                y_label = 'Star Mass',
                y_values = y_values,
                label = plot_label)

#StarCount graph
x_values = rings.keys() #the ring indexes e.g. 0, 1, 2, 3...
y_values = [count for (count, mass, density) in rings.values()]
plot_star_rings(filename = 'star_count_plot.pdf',
                title = 'Scatter Plot of Star Radius (in rings) vs. Star Count',
                x_label = 'Rings',
                x_values = x_values,
                y_label = 'Star Count',
                y_values = y_values,
                label = plot_label)  

#accumlated star mass graph plot
ring_indexes = rings.keys() #the ring indexes e.g. 0, 1, 2, 3...
ring_star_mass = [mass for (count, mass, density) in rings.values()]
accumlated_star_mass = []
for i in range(len(ring_star_mass)):
    if i == 0:
        accumlated_star_mass.append(ring_star_mass[i])
    else:
        accumlated_star_mass.append(ring_star_mass[i]+accumlated_star_mass[i-1])
plot_star_rings(filename = 'accumulated_star_mass_plot.pdf',
                title = 'Scatter Plot of Star Radius (in rings) vs. Accumulated Star Mass',
                x_label = 'Rings',
                x_values = ring_indexes,
                y_label = 'Accumulated Star Mass',
                y_values = accumlated_star_mass,
                label = plot_label)   
#Star_density plot
lot_label = 'Ring Thickness = {}, \n Max Radius = {} KPC'.format(slice_thickness, max_radius)
x_values = rings.keys() #the ring indexes e.g. 0, 1, 2, 3...
y_values = [density for (count, mass, density) in rings.values()]
plot_star_rings(filename = 'star_mass_density_plot.pdf',
                title = 'Scatter Plot of Star Radius (in rings) vs. Star Density',
                x_label = 'Rings',
                x_values = x_values,
                y_label = 'Star Density',
                y_values = y_values,
                label = plot_label)