

import p_classes.class_particles as class_particle
print("Start")

#...Dashboard: (choose from dm, stars, and gas)
component = 'dm'

h_obj = class_particle.component_graphs()

#radius ranges from ~ 0 - 17.7
h_obj.mask(radius_inner = 0, radius_outer = 17.7, component = component)

h_obj.data_arrange(component=component,var_1_name = 'x', var_2_name = 'y', color = 'blue', symbol = '^', size = 50, weight = 'mass')

h_obj.total_mass()

h_obj.particles_ring(component=component)



print("End")





