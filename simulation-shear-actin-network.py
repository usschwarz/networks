#Script to set up and shear a cytoskeletal network with actin filaments
#used to generate data for Figure 4H
#author: Julia Jaeger

import numpy as np
import matplotlib.pyplot as plt
import readdy
import os

#function to extract the stress from readdy trajectory data

def extract_quantities_own(traj,mean_xy,mean_yx):
    
	t, forces = traj.read_observable_forces()
	t, types, ids, pos = traj.read_observable_particles()

	tau_xy=np.zeros(len(pos))
	tau_yx=np.zeros(len(pos))

	f_tot_x=np.zeros(len(pos))
	f_tot_y=np.zeros(len(pos))

	for time in range(len(pos)):
		for particle in range(len(pos[0])):
			f_i=forces[time][particle]
			r_i=pos[time][particle]

			tau_xy[time]+=f_i[1]*r_i[0]
			tau_yx[time]+=f_i[0]*r_i[1]

			f_tot_x[time]+=f_i[0]
			f_tot_y[time]+=f_i[1]

	equil=int(len(tau_xy)*0.2)   #check in plot how many to use (equilibration)
	mean_xy.append(np.mean(tau_xy[equil:]))
	mean_yx.append(np.mean(tau_yx[equil:]))

	return f_tot_x,f_tot_y

#function to extract the angles of the actin filaments

def angles(p1,p2):
	dx=np.absolute(p1[0]-p2[0])
	dy=np.absolute(p1[1]-p2[1])
	dr=np.sqrt(dx*dx+dy*dy)
	dz=np.absolute(p1[2]-p2[2])
	z_angle=np.arctan(dz/dr)
	plane_angle=np.arctan(dy/dx)
	return z_angle,plane_angle


def extract_angles(traj,angle_z,angle_plane,actin_positions):

	t, types, ids, pos = traj.read_observable_particles()

	type_names=[]
	for t_id in types[-1]:
	    	type_names.append(traj.species_name(t_id))
	positions=pos[-1]

	ind=0
	for i in range(46):

		while type_names[ind]!='pointedcap' and type_names[ind]!='pointed':
			ind+=1

		if type_names[ind]=='pointedcap':
            		pc=True
		elif type_names[ind]=='pointed':
            		pc=False
		s=ind
		ind+=1
		while type_names[ind]=='core':
            		ind+=1
		if type_names[ind]=='barbedcap':
            		bc=True
		elif type_names[ind]=='barbed':
            		bc=False
		e=ind

		actin_positions.append([np.mean(positions[s:e].T[0]),np.mean(positions[s:e].T[1]),np.mean(positions[s:e].T[2])])
		z_a,plane_a=angles(positions[s],positions[e])
		angle_z.append(z_a)
		angle_plane.append(plane_a)

		ind+=1

#definition of dissociation rates for readdy

def pointed_rate(topology):
	types = [topology.particle_type_of_vertex(v) for v in topology.get_graph().get_vertices()]
	n = len(topology.get_graph().get_vertices())
	if n>2 and ("pointed" in types):
		r=0.0000008
	else:
		r=0.0 
	return r

def pointed_dissociation(topology):
	recipe = readdy.StructuralReactionRecipe(topology)
	    
	vertex_to_separate = None
	for index, vertex in enumerate(topology.get_graph().get_vertices()):
		if "pointed" == topology.particle_type_of_vertex(vertex):
			vertex_to_separate = index
			v=vertex
	assert vertex_to_separate is not None
	      
	# separate particle
	vertex_neighbour= None
	n_neighbors=0
	for neighbor in v:
		n_neighbors += 1
		vertex_neighbour = neighbor.get().particle_index
	assert vertex_neighbour is not None

	# change adjacent type to head
	recipe.change_particle_type(vertex_to_separate, "Gactin") 
	recipe.change_particle_type(vertex_neighbour, "pointed")
	recipe.separate_vertex(vertex_to_separate)
	return recipe

def barbed_rate(topology):
	types = [topology.particle_type_of_vertex(v) for v in topology.get_graph().get_vertices()]
	n = len(topology.get_graph().get_vertices())
	if n>2 and ("barbed" in types):
		r=0.0000014
	else:
		r=0.
	return r

def barbed_dissociation(topology):
	recipe = readdy.StructuralReactionRecipe(topology)
	vertex_to_separate = None
	for index, vertex in enumerate(topology.get_graph().get_vertices()):
		if "barbed" == topology.particle_type_of_vertex(vertex):
	    		vertex_to_separate = index
	    		v=vertex
	assert vertex_to_separate is not None
    
	# separate particle

	n_neighbors=0
	for neighbor in v:
		n_neighbors += 1
		vertex_neighbour = neighbor.get().particle_index
	assert vertex_neighbour is not None
    
	# change adjacent type to head
	recipe.change_particle_type(vertex_neighbour, "barbed")
	recipe.change_particle_type(vertex_to_separate, "Gactin") 
	recipe.separate_vertex(vertex_to_separate)
	return recipe

def pointedcap_rate(topology):
	types = [topology.particle_type_of_vertex(v) for v in topology.get_graph().get_vertices()]
	n = len(topology.get_graph().get_vertices())
	if n>2 and ("pointedcap" in types):
		r=0.0000008
	else:
		r=0.
	return r
    
def pointedcap_dissociation(topology):
	recipe = readdy.StructuralReactionRecipe(topology)
	vertex_to_separate = None
	for index, vertex in enumerate(topology.get_graph().get_vertices()):
		if "pointedcap" == topology.particle_type_of_vertex(vertex):
	    		vertex_to_separate = index
	    		v=vertex
	assert vertex_to_separate is not None
    
	# separate particle
	n_neighbors=0
	for neighbor in v:
	    	n_neighbors += 1
	    	vertex_neighbour = neighbor.get().particle_index
	assert vertex_neighbour is not None
    
	# change adjacent type to head
	recipe.change_particle_type(vertex_neighbour, "pointed")
	recipe.change_particle_type(vertex_to_separate, "tropomodulin")
	recipe.separate_vertex(vertex_to_separate)
	return recipe

def barbedcap_rate(topology):
	types = [topology.particle_type_of_vertex(v) for v in topology.get_graph().get_vertices()]
	n = len(topology.get_graph().get_vertices())
	if n>2 and ("barbedcap" in types):
		r=0.0000014
	else:
		r=0.
	return r
    
def barbedcap_dissociation(topology):
	recipe = readdy.StructuralReactionRecipe(topology)
	vertex_to_separate = None
	for index, vertex in enumerate(topology.get_graph().get_vertices()):
		if "barbedcap" == topology.particle_type_of_vertex(vertex):
	    		vertex_to_separate = index
	    		v=vertex
	assert vertex_to_separate is not None
    
	# separate particle
	n_neighbors=0
	for neighbor in v:
	    	n_neighbors += 1
	    	vertex_neighbour = neighbor.get().particle_index
	assert vertex_neighbour is not None
    
	# change adjacent type to head
	recipe.change_particle_type(vertex_neighbour, "barbed")
	recipe.change_particle_type(vertex_to_separate, "adducin")
	recipe.separate_vertex(vertex_to_separate)
	return recipe

def random_sample(N, z_min, z_max):
	positions=np.zeros((N,3))
	for i in range(N):
		positions[i][0]=np.random.rand()*2*x_base - x_base
		positions[i][1]=np.random.rand()*y_base - y_base* 0.5
		positions[i][2]=z_min+np.random.rand()*(z_max-z_min)
	return positions


#function to set up and simulte the initial equilibration of the cytoskeleton

def set_up_and_run_first(x_base,y_base,out_file,cps,writeout,checkpoints,runtime):

	system = readdy.ReactionDiffusionSystem(box_size=[2*x_base*1.2, y_base*1.2, 100+100])
	system.periodic_boundary_conditions = [False, False, False]
	volume=0.001*2*x_base*1.1*0.001*y_base*1.1*0.1 #um^3  (constricted by potential in z-direction)
	system.temperature=300

	z_pos=-40

	g_actins = 20#126
	kahrps=0
	tropomodulins= 5
	adducins= 5


	# In[4]:

	epsilon_k_k=system.kbt*0.1

	epsilon_k_sp=epsilon_k_k
	epsilon_k_a=epsilon_k_k


	#definition of particles with diffusion constants

	system.add_species("Gactin", 71.5 * readdy.units.um**2 / readdy.units.s)
	system.add_species("adducin", 51.07 * readdy.units.um**2 / readdy.units.s)
	system.add_species("tropomodulin", 29.59 * readdy.units.um**2 / readdy.units.s)
	system.add_species("kahrp", 76.6 * readdy.units.um**2 / readdy.units.s) 

	system.topologies.add_type("filament") #reduced diffusion
	system.add_topology_species("core", 1 * readdy.units.um**2 / readdy.units.s)  
	system.add_topology_species("pointed", 1 * readdy.units.um**2 / readdy.units.s)   
	system.add_topology_species("barbed", 0.53 * readdy.units.um**2 / readdy.units.s)  
	system.add_topology_species("pointedcap", 1 * readdy.units.um**2 / readdy.units.s) 
	system.add_topology_species("barbedcap", 0.53 * readdy.units.um**2 / readdy.units.s)  

	system.topologies.add_type("spectrin")
	system.add_topology_species("spectrincore", 41.1 * readdy.units.um**2 / readdy.units.s)
	system.add_topology_species("spectrinedge", 41.1 * readdy.units.um**2 / readdy.units.s)
	system.add_topology_species("spectrinmidpoint", 0.53 * readdy.units.um**2 / readdy.units.s)
	system.add_topology_species("spectrinbinding", 41.1 * readdy.units.um**2 / readdy.units.s)

	system.add_species("anchor", 0.53 * readdy.units.um**2 / readdy.units.s)  

	# In[5]:


	#add box and plane potential

	for i in ["Gactin","adducin","kahrp","tropomodulin","core","pointed","barbed","pointedcap","barbedcap","spectrincore","spectrinbinding","spectrinedge","spectrinmidpoint","anchor"]: 
	    system.potentials.add_box(
		particle_type=i, force_constant=100, origin=[-x_base*1.1, -y_base*1.1*0.5, -50], 
		extent=[2*x_base*1.1, y_base*1.1,100]
	    )
	for i in ["core","pointed","barbed","pointedcap","barbedcap","spectrinmidpoint"]:   #"spectrinedge"
	    system.potentials.add_box(
		particle_type=i, force_constant=10, origin=[-x_base*1.1, -y_base*1.1*0.5, -42], #confinement in plane in middle of box
		extent=[2*x_base*1.1, y_base*1.1, 10]
	    )


	# In[6]:


	particle_radii = {"Gactin": 3., "adducin": 4.2, "kahrp": 2.8, "tropomodulin": 7.25, "core": 3., "pointed": 3.,"pointedcap": 7.25, "barbed": 3., "barbedcap": 4.2, "spectrincore": 5.26/2.0, "spectrinedge": 5.26/2.0,"spectrinmidpoint":5.26/2.0,"spectrinbinding":5.26/2.0,"anchor": 5.0} #radius in nm

	epsilon_const=14*4.184

	sigma_sp_a=(particle_radii["spectrincore"]+particle_radii["Gactin"])#np.power(2,5/6.)*
	sigma_sp=(particle_radii["spectrincore"]+particle_radii["spectrincore"])
	sigma_a=(particle_radii["Gactin"]+particle_radii["Gactin"])
	sigma_k_a=(particle_radii["kahrp"]+particle_radii["Gactin"])
	sigma_k_sp=(particle_radii["kahrp"]+particle_radii["spectrincore"])
	sigma_k_k=(particle_radii["kahrp"]+particle_radii["kahrp"])

	bond_force_spectrin=36*np.power(2,2./3)*epsilon_const/np.power(sigma_sp,2)
	angle_force_spectrin=2.62
	repulsion_force_spectrin=36*np.power(2,2./3)*epsilon_const/np.power(sigma_sp,2)

	bond_force=36*np.power(2,2./3)*epsilon_const/np.power(sigma_a,2)
	angle_force=4280.
	repulsion_force=36*np.power(2,2./3)*epsilon_const/np.power(sigma_a,2)

	repulsion_force_sp_a=36*np.power(2,2./3)*epsilon_const/np.power(sigma_sp_a,2)
	repulsion_force_k_a=36*np.power(2,2./3)*epsilon_const/np.power(sigma_k_a,2)
	repulsion_force_k_sp=36*np.power(2,2./3)*epsilon_const/np.power(sigma_k_sp,2)


	#set up inter particle potentials
	#direction does not matter (AB=BA)
	system.topologies.configure_harmonic_bond("core", "core", force_constant=bond_force, length=particle_radii["core"]+particle_radii["core"])
	system.topologies.configure_harmonic_bond("core", "pointed", force_constant=bond_force, length=particle_radii["core"]+particle_radii["pointed"])
	system.topologies.configure_harmonic_bond("core", "barbed", force_constant=bond_force, length=particle_radii["core"]+particle_radii["barbed"])
	system.topologies.configure_harmonic_bond("core", "pointedcap", force_constant=bond_force, length=particle_radii["core"]+particle_radii["pointedcap"])
	system.topologies.configure_harmonic_bond("core", "barbedcap", force_constant=bond_force, length=particle_radii["core"]+particle_radii["barbedcap"])
	system.topologies.configure_harmonic_bond("pointed", "barbed", force_constant=bond_force, length=particle_radii["pointed"]+particle_radii["barbed"])
	system.topologies.configure_harmonic_bond("pointedcap", "barbed", force_constant=bond_force, length=particle_radii["pointedcap"]+particle_radii["barbed"])
	system.topologies.configure_harmonic_bond("pointed", "barbedcap", force_constant=bond_force, length=particle_radii["pointed"]+particle_radii["barbedcap"])

	system.topologies.configure_harmonic_angle("core", "core", "core", force_constant=angle_force, equilibrium_angle=np.pi)
	system.topologies.configure_harmonic_angle("pointed", "core", "core", force_constant=angle_force, equilibrium_angle=np.pi)
	system.topologies.configure_harmonic_angle("barbed", "core", "core", force_constant=angle_force, equilibrium_angle=np.pi)
	system.topologies.configure_harmonic_angle("pointedcap", "core", "core", force_constant=angle_force, equilibrium_angle=np.pi)
	system.topologies.configure_harmonic_angle("barbedcap", "core", "core", force_constant=angle_force, equilibrium_angle=np.pi)
	system.topologies.configure_harmonic_angle("pointed", "core", "barbed", force_constant=angle_force, equilibrium_angle=np.pi)
	system.topologies.configure_harmonic_angle("barbed", "core", "pointed", force_constant=angle_force, equilibrium_angle=np.pi)
	system.topologies.configure_harmonic_angle("pointedcap", "core", "barbed", force_constant=angle_force, equilibrium_angle=np.pi)
	system.topologies.configure_harmonic_angle("barbedcap", "core", "pointed", force_constant=angle_force, equilibrium_angle=np.pi)

	system.topologies.configure_harmonic_bond("spectrincore", "spectrincore", force_constant=bond_force_spectrin, length=particle_radii["spectrincore"]+particle_radii["spectrincore"])
	system.topologies.configure_harmonic_bond("spectrincore", "spectrinedge", force_constant=bond_force_spectrin, length=particle_radii["spectrincore"]+particle_radii["spectrinedge"])
	system.topologies.configure_harmonic_bond("spectrincore", "spectrinbinding", force_constant=bond_force_spectrin, length=particle_radii["spectrincore"]+particle_radii["spectrinbinding"])
	system.topologies.configure_harmonic_bond("spectrincore", "spectrinmidpoint", force_constant=bond_force_spectrin, length=particle_radii["spectrincore"]+particle_radii["spectrinmidpoint"])
	system.topologies.configure_harmonic_bond("spectrinbinding", "spectrinbinding", force_constant=bond_force_spectrin, length=particle_radii["spectrinbinding"]+particle_radii["spectrinbinding"])
	system.topologies.configure_harmonic_bond("spectrinmidpoint", "spectrinbinding", force_constant=bond_force_spectrin, length=particle_radii["spectrinmidpoint"]+particle_radii["spectrinbinding"])
	system.topologies.configure_harmonic_bond("spectrinmidpoint", "spectrinmidpoint", force_constant=bond_force_spectrin, length=particle_radii["spectrinmidpoint"]+particle_radii["spectrinmidpoint"])

	system.topologies.configure_harmonic_angle("spectrincore", "spectrincore", "spectrincore", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system.topologies.configure_harmonic_angle("spectrinedge", "spectrincore", "spectrincore", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system.topologies.configure_harmonic_angle("spectrincore", "spectrinbinding", "spectrincore", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system.topologies.configure_harmonic_angle("spectrincore", "spectrincore", "spectrinbinding", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system.topologies.configure_harmonic_angle("spectrincore", "spectrinbinding", "spectrinbinding", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system.topologies.configure_harmonic_angle("spectrinbinding", "spectrinbinding", "spectrinbinding", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system.topologies.configure_harmonic_angle("spectrincore", "spectrincore", "spectrinmidpoint", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system.topologies.configure_harmonic_angle("spectrincore", "spectrinmidpoint", "spectrincore", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system.topologies.configure_harmonic_angle("spectrincore", "spectrinmidpoint", "spectrinmidpoint", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system.topologies.configure_harmonic_angle("spectrinbinding", "spectrinmidpoint", "spectrinmidpoint", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system.topologies.configure_harmonic_angle("spectrinbinding", "spectrinbinding", "spectrinmidpoint", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system.topologies.configure_harmonic_angle("spectrinbinding", "spectrincore", "spectrinbinding", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system.topologies.configure_harmonic_angle("spectrinbinding", "spectrincore", "spectrinmidpoint", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system.topologies.configure_harmonic_angle("spectrincore", "spectrinbinding", "spectrincore", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)

	#harmonic repulsion of freely diffusing particles
	all_pairs = [("Gactin","Gactin"), ("adducin","adducin"), ("tropomodulin", "tropomodulin"), ("Gactin", "adducin"), ("Gactin", "tropomodulin"), ("adducin", "tropomodulin")]
	for pair in all_pairs:
	    distance = particle_radii[pair[0]] + particle_radii[pair[1]]
	    system.potentials.add_harmonic_repulsion(pair[0], pair[1], force_constant=repulsion_force, 
		                                     interaction_distance=distance)

	#KAHRP with actin
	all_pairs = [("kahrp","Gactin"), ("kahrp", "adducin"), ("kahrp", "tropomodulin"), ("pointed", "kahrp"), ("pointedcap", "kahrp"), ("barbed", "kahrp"), ("barbedcap", "kahrp")]
	for pair in all_pairs:
	    distance = particle_radii[pair[0]] + particle_radii[pair[1]]
	    system.potentials.add_harmonic_repulsion(pair[0], pair[1], force_constant=repulsion_force_k_a, 
		                                     interaction_distance=distance)

	#KAHRP with spectrin
	all_pairs = [("kahrp","spectrincore"),("kahrp","spectrinedge"), ("spectrinmidpoint","kahrp")] #spectrinbinding -> LJ
	for pair in all_pairs:
	    distance = particle_radii[pair[0]] + particle_radii[pair[1]]
	    system.potentials.add_harmonic_repulsion(pair[0], pair[1], force_constant=repulsion_force_k_a, 
		                                     interaction_distance=distance)
	    
	#harmonic repulsion of freely diffusing particles with core
	all_pairs = [("core","core"),("Gactin","core"),("Gactin","pointedcap"),("Gactin","barbedcap"), ("adducin","core"), ("adducin","pointed"), ("adducin","pointedcap"), ("tropomodulin", "core"), ("tropomodulin", "barbed"), ("tropomodulin", "barbedcap"), ("tropomodulin", "pointedcap"), ("adducin", "barbedcap")]
	for pair in all_pairs:
	    distance = particle_radii[pair[0]] + particle_radii[pair[1]]
	    system.potentials.add_harmonic_repulsion(pair[0], pair[1], force_constant=repulsion_force, 
		                                     interaction_distance=distance)

	#harmonic repulsion of actin filaments
	all_pairs = [("core","core"),("pointed","core"),("pointedcap","core"),("barbedcap","core"),("barbed","core"),("barbed","pointed"),("barbedcap","pointed"),("pointedcap","pointed"),("pointed","pointed"),("barbed","pointedcap"),("barbed","barbed"),("barbed","barbedcap"),("pointedcap","barbedcap"),("pointedcap","pointedcap"),("barbedcap","barbedcap")]
	for pair in all_pairs:
	    distance = particle_radii[pair[0]] + particle_radii[pair[1]]
	    system.potentials.add_harmonic_repulsion(pair[0], pair[1], force_constant=repulsion_force, 
		                                     interaction_distance=distance)

	#harmonic repulsion of freely diffusing particles with spectrin
	all_pairs = [("Gactin","spectrincore"), ("adducin","spectrincore"), ("tropomodulin", "spectrincore"),("Gactin","spectrinbinding"), ("adducin","spectrinbinding"), ("tropomodulin", "spectrinbinding"),("Gactin","spectrinedge"), ("adducin","spectrinedge"), ("tropomodulin", "spectrinedge"),("Gactin","spectrinmidpoint"), ("adducin","spectrinmidpoint"), ("tropomodulin", "spectrinmidpoint")]
	for pair in all_pairs:
	    distance = particle_radii[pair[0]] + particle_radii[pair[1]]
	    system.potentials.add_harmonic_repulsion(pair[0], pair[1], force_constant=repulsion_force_sp_a, 
		                                     interaction_distance=distance)
	    
	#harmonic repulsion of spectrin core 
	all_pairs = [("core","spectrincore"), ("pointed", "spectrincore"), ("pointedcap", "spectrincore"), ("barbed", "spectrincore"), ("barbedcap", "spectrincore")]
	for pair in all_pairs:
	    distance = particle_radii[pair[0]] + particle_radii[pair[1]]
	    system.potentials.add_harmonic_repulsion(pair[0], pair[1], force_constant=repulsion_force_sp_a, 
		                                     interaction_distance=distance)

	#harmonic repulsion of spectrin heads
	all_pairs = [("spectrinedge","spectrinedge")]
	for pair in all_pairs:
	    distance = (particle_radii[pair[0]] + particle_radii[pair[1]])
	    system.potentials.add_harmonic_repulsion(pair[0], pair[1], force_constant=repulsion_force_spectrin*2, 
		                                     interaction_distance=distance*2)

	#harmonic repulsion of spectrin midpoints
	all_pairs = [("spectrinmidpoint","core"),("spectrinmidpoint","barbed"),("spectrinmidpoint","barbedcap"),("spectrinmidpoint","pointed"),("spectrinmidpoint","pointedcap")]
	for pair in all_pairs:
	    distance = (particle_radii[pair[0]] + particle_radii[pair[1]])
	    system.potentials.add_harmonic_repulsion(pair[0], pair[1], force_constant=repulsion_force_sp_a, 
		                                     interaction_distance=distance)

	all_pairs = [("spectrinbinding","core"),("spectrinbinding","barbed"),("spectrinbinding","barbedcap"), ("spectrinedge","barbed"), ("spectrinedge","barbedcap"),("spectrinbinding","pointed"),("spectrinedge","pointed"),("spectrinedge","pointedcap"),("spectrinbinding","pointedcap")]
	for pair in all_pairs:
	    distance = (particle_radii[pair[0]] + particle_radii[pair[1]])
	    system.potentials.add_harmonic_repulsion(pair[0], pair[1], force_constant=repulsion_force_sp_a, 
		                                     interaction_distance=distance)

	all_pairs = [("spectrinbinding","spectrinbinding"),("spectrinbinding","spectrincore"),("spectrincore","spectrincore"),("spectrincore","spectrinedge"), ("spectrinmidpoint","spectrinmidpoint"), ("spectrinmidpoint","spectrinedge"), ("spectrinbinding","spectrinedge"), ("spectrinmidpoint","spectrincore"),("spectrinmidpoint","spectrinbinding")]
	for pair in all_pairs:
	    distance = (particle_radii[pair[0]] + particle_radii[pair[1]])
	    system.potentials.add_harmonic_repulsion(pair[0], pair[1], force_constant=repulsion_force_spectrin, 
		                                     interaction_distance=distance)


	#Interaction of Spectrinedge with actin (core, pointed, barbed)

	system.potentials.add_lennard_jones("spectrinedge","core", m=12, n=6, cutoff=20, shift=True, epsilon=epsilon_const, sigma=sigma_sp_a)
	system.potentials.add_lennard_jones("spectrinedge","anchor", m=12, n=6, cutoff=20, shift=True, epsilon=epsilon_const*5, sigma=np.power(2,5/6.)*(particle_radii["spectrinedge"]+particle_radii["anchor"])/2.)
	
	#Binding to filament

	binding_rate_pointed=0.0000686  #in 1/ns  #0.0000068 
	binding_rate_barbed=0.000597 #in 1/ns  #0.00016
	binding_rate_tropomodulin=24324  #in 1/ns
	binding_rate_adducin=0.00329 #in 1/ns
	binding_radius=6.  #to middle of G-actin

	system.topologies.add_spatial_reaction(
	    "BindPointed: filament(pointed) + (Gactin) -> filament(core--pointed)", 
	    rate=binding_rate_pointed, radius=binding_radius
	)
	system.topologies.add_spatial_reaction(
	    "BindBarbed: filament(barbed) + (Gactin) -> filament(core--barbed)", 
	    rate=binding_rate_barbed, radius=binding_radius
	)
	system.topologies.add_spatial_reaction(
	    "BindPointedCap: filament(pointed) + (tropomodulin) -> filament(core--pointedcap)", 
	    rate=binding_rate_tropomodulin, radius=binding_radius
	)
	system.topologies.add_spatial_reaction(
	    "BindBarbedCap: filament(barbed) + (adducin) -> filament(core--barbedcap)", 
	    rate=binding_rate_adducin, radius=binding_radius
	)


	#dissociation of end beads

	system.topologies.add_structural_reaction("pointed_dissociation","filament", pointed_dissociation, pointed_rate)
	system.topologies.add_structural_reaction("barbed_dissociation", "filament", barbed_dissociation, barbed_rate)
	system.topologies.add_structural_reaction("pointedcap_dissociation", "filament", pointedcap_dissociation, pointedcap_rate)
	system.topologies.add_structural_reaction("barbedcap_dissociation", "filament", barbedcap_dissociation, barbedcap_rate)

	simulation = system.simulation(kernel="CPU")


	#set up the cytoskeletal filaments

	left_border=-x_base
	x_coords=[left_border,left_border+0.5*a, left_border+a, left_border+1.5*a, left_border+2.0*a, left_border+2.5*a, left_border+3.0*a,left_border+3.5*a, left_border+4.0*a, left_border+4.5*a, left_border+5.0*a, left_border+5.5*a, left_border+6.0*a, left_border+6.5*a, left_border+7.0*a]
	y_coords=[-y_base*0.5,-y_base*0.5+h,-y_base*0.5+2*h,-y_base*0.5+3*h,-y_base*0.5+4*h,-y_base*0.5+5*h,-y_base*0.5+6*h,-y_base*0.5+7*h,y_base*0.5]

	actin_angles=np.random.rand(46)*2*np.pi

	#initialize actin filaments (no periodicity needed here)
	def periodic_x(x):
		return x
		#if x>=x_coords[0] and x<=x_coords[4]:
	#		x_new=x
	#	elif x>x_coords[4]:
	#		x_new=x_coords[0]+x-x_coords[4]
	#	elif x<x_coords[4]:
	#		x_new=x_coords[4]+(x-x_coords[0])
	#	return x_new

	def periodic_y(y):
		return y
		#if y>=y_coords[0] and y<=y_coords[4]:
	#		y_new=y
	#	elif y>y_coords[4]:
	#		y_new=y_coords[0]+y-y_coords[4]
	#	elif y<y_coords[4]:
	#		y_new=y_coords[4]+(y-y_coords[0])
	#	return y_new

	def initialize_actin(x,y,z,alpha):
		init_top_pos = np.array([
		    [periodic_x(x+(-25.)*np.cos(alpha)),periodic_y(y+(-25.)*np.sin(alpha)),z],
		    [periodic_x(x+(-15.)*np.cos(alpha)),periodic_y(y+(-15.)*np.sin(alpha)),z],
		    [periodic_x(x+(-9.)*np.cos(alpha)),periodic_y(y+(-9.)*np.sin(alpha)),z],
		    [periodic_x(x+(-3.)*np.cos(alpha)),periodic_y(y+(-3.)*np.sin(alpha)),z],
		    [periodic_x(x+(3.)*np.cos(alpha)),periodic_y(y+(3.)*np.sin(alpha)),z],
		    [periodic_x(x+(9.)*np.cos(alpha)),periodic_y(y+(9.)*np.sin(alpha)),z],
		    [periodic_x(x+(15.)*np.cos(alpha)),periodic_y(y+(15.)*np.sin(alpha)),z],
		    [periodic_x(x+(21.)*np.cos(alpha)),periodic_y(y+(21.)*np.sin(alpha)),z]
		])
		top = simulation.add_topology("filament", ["pointedcap", "core", "core", "core", "core", "core", "core", "barbedcap"], init_top_pos)
		top.get_graph().add_edge(0, 1)
		top.get_graph().add_edge(1, 2)
		top.get_graph().add_edge(2, 3)
		top.get_graph().add_edge(3, 4)
		top.get_graph().add_edge(4, 5)	
		top.get_graph().add_edge(5, 6)	
		top.get_graph().add_edge(6, 7)	


	#46 actins

	initialize_actin(x_coords[1],y_coords[1],z_pos,actin_angles[0])
	initialize_actin(x_coords[3],y_coords[1],z_pos,actin_angles[1])
	initialize_actin(x_coords[5],y_coords[1],z_pos,actin_angles[2])
	initialize_actin(x_coords[7],y_coords[1],z_pos,actin_angles[3])
	initialize_actin(x_coords[9],y_coords[1],z_pos,actin_angles[4])
	initialize_actin(x_coords[11],y_coords[1],z_pos,actin_angles[5])
	initialize_actin(x_coords[13],y_coords[1],z_pos,actin_angles[6])

	initialize_actin(x_coords[2],y_coords[2],z_pos,actin_angles[7])
	initialize_actin(x_coords[4],y_coords[2],z_pos,actin_angles[8])
	initialize_actin(x_coords[6],y_coords[2],z_pos,actin_angles[9])
	initialize_actin(x_coords[8],y_coords[2],z_pos,actin_angles[10])
	initialize_actin(x_coords[10],y_coords[2],z_pos,actin_angles[11])
	initialize_actin(x_coords[12],y_coords[2],z_pos,actin_angles[12])

	initialize_actin(x_coords[1],y_coords[3],z_pos,actin_angles[13])
	initialize_actin(x_coords[3],y_coords[3],z_pos,actin_angles[14])
	initialize_actin(x_coords[5],y_coords[3],z_pos,actin_angles[15])
	initialize_actin(x_coords[7],y_coords[3],z_pos,actin_angles[16])
	initialize_actin(x_coords[9],y_coords[3],z_pos,actin_angles[17])
	initialize_actin(x_coords[11],y_coords[3],z_pos,actin_angles[18])
	initialize_actin(x_coords[13],y_coords[3],z_pos,actin_angles[19])

	initialize_actin(x_coords[2],y_coords[4],z_pos,actin_angles[20])
	initialize_actin(x_coords[4],y_coords[4],z_pos,actin_angles[21])
	initialize_actin(x_coords[6],y_coords[4],z_pos,actin_angles[22])
	initialize_actin(x_coords[8],y_coords[4],z_pos,actin_angles[23])
	initialize_actin(x_coords[10],y_coords[4],z_pos,actin_angles[24])
	initialize_actin(x_coords[12],y_coords[4],z_pos,actin_angles[25])	

	initialize_actin(x_coords[1],y_coords[5],z_pos,actin_angles[26])
	initialize_actin(x_coords[3],y_coords[5],z_pos,actin_angles[27])
	initialize_actin(x_coords[5],y_coords[5],z_pos,actin_angles[28])
	initialize_actin(x_coords[7],y_coords[5],z_pos,actin_angles[29])
	initialize_actin(x_coords[9],y_coords[5],z_pos,actin_angles[30])
	initialize_actin(x_coords[11],y_coords[5],z_pos,actin_angles[31])
	initialize_actin(x_coords[13],y_coords[5],z_pos,actin_angles[32])

	initialize_actin(x_coords[2],y_coords[6],z_pos,actin_angles[33])
	initialize_actin(x_coords[4],y_coords[6],z_pos,actin_angles[34])
	initialize_actin(x_coords[6],y_coords[6],z_pos,actin_angles[35])
	initialize_actin(x_coords[8],y_coords[6],z_pos,actin_angles[36])
	initialize_actin(x_coords[10],y_coords[6],z_pos,actin_angles[37])
	initialize_actin(x_coords[12],y_coords[6],z_pos,actin_angles[38])

	initialize_actin(x_coords[1],y_coords[7],z_pos,actin_angles[39])
	initialize_actin(x_coords[3],y_coords[7],z_pos,actin_angles[40])
	initialize_actin(x_coords[5],y_coords[7],z_pos,actin_angles[41])
	initialize_actin(x_coords[7],y_coords[7],z_pos,actin_angles[42])
	initialize_actin(x_coords[9],y_coords[7],z_pos,actin_angles[43])
	initialize_actin(x_coords[11],y_coords[7],z_pos,actin_angles[44])
	initialize_actin(x_coords[13],y_coords[7],z_pos,actin_angles[45])




	def spectrin_initialization_periodic(actin_x_start,actin_y_start,actin_x_end,actin_y_end,pos_start,pos_end,z_plane,actin_angle_start,actin_angle_end,choices):
		d_sp=7.
		
		#move start and end according to angle
		if actin_angle_start<=np.pi/4. or 7*np.pi/4.<actin_angle_start:
			if pos_start-4<0.001:
				startx=actin_x_start+(15.)*np.cos(actin_angle_start)-(d_sp)*np.sin(actin_angle_start)
				starty=actin_y_start+(15.)*np.sin(actin_angle_start)+(d_sp)*np.cos(actin_angle_start)
			elif pos_start-5<0.001:
				startx=actin_x_start+(15.)*np.cos(actin_angle_start)-(-d_sp)*np.sin(actin_angle_start)
				starty=actin_y_start+(15.)*np.sin(actin_angle_start)+(-d_sp)*np.cos(actin_angle_start)
			elif pos_start-6<0.001:
				startx=actin_x_start+(0.)*np.cos(actin_angle_start)-(-d_sp)*np.sin(actin_angle_start)
				starty=actin_y_start+(0.)*np.sin(actin_angle_start)+(-d_sp)*np.cos(actin_angle_start)
		elif (np.pi/4.<actin_angle_start and actin_angle_start<=3*np.pi/4.): 
			#special case of anchor			
			if pos_start-2<0.001:
				startx=actin_x_start+(0.)*np.cos(actin_angle_start)-(d_sp)*np.sin(actin_angle_start)
				starty=actin_y_start+(0.)*np.sin(actin_angle_start)+(d_sp)*np.cos(actin_angle_start)
			elif pos_start-4<0.001:
				startx=actin_x_start+(15.)*np.cos(actin_angle_start)-(-d_sp)*np.sin(actin_angle_start)
				starty=actin_y_start+(15.)*np.sin(actin_angle_start)+(-d_sp)*np.cos(actin_angle_start)
			elif pos_start-5<0.001:
				startx=actin_x_start+(0.)*np.cos(actin_angle_start)-(-d_sp)*np.sin(actin_angle_start)
				starty=actin_y_start+(0.)*np.sin(actin_angle_start)+(-d_sp)*np.cos(actin_angle_start)
			elif pos_start-6<0.001:
				startx=actin_x_start+(-15.)*np.cos(actin_angle_start)-(-d_sp)*np.sin(actin_angle_start)
				starty=actin_y_start+(-15.)*np.sin(actin_angle_start)+(-d_sp)*np.cos(actin_angle_start)
		elif (3*np.pi/4.<actin_angle_start and actin_angle_start<=5*np.pi/4.): 
			if pos_start-4<0.001:
				startx=actin_x_start+(0.)*np.cos(actin_angle_start)-(-d_sp)*np.sin(actin_angle_start)
				starty=actin_y_start+(0.)*np.sin(actin_angle_start)+(-d_sp)*np.cos(actin_angle_start)
			elif pos_start-5<0.001:
				startx=actin_x_start+(-15.)*np.cos(actin_angle_start)-(-d_sp)*np.sin(actin_angle_start)
				starty=actin_y_start+(-15.)*np.sin(actin_angle_start)+(-d_sp)*np.cos(actin_angle_start)
			elif pos_start-6<0.001:
				startx=actin_x_start+(-15.)*np.cos(actin_angle_start)-(d_sp)*np.sin(actin_angle_start)
				starty=actin_y_start+(-15.)*np.sin(actin_angle_start)+(d_sp)*np.cos(actin_angle_start)	
		elif (5*np.pi/4.<actin_angle_start and actin_angle_start<=7*np.pi/4.): 
			if pos_start-4<0.001:
				startx=actin_x_start+(-15.)*np.cos(actin_angle_start)-(d_sp)*np.sin(actin_angle_start)
				starty=actin_y_start+(-15.)*np.sin(actin_angle_start)+(d_sp)*np.cos(actin_angle_start)
			elif pos_start-5<0.001:
				startx=actin_x_start+(0.)*np.cos(actin_angle_start)-(d_sp)*np.sin(actin_angle_start)
				starty=actin_y_start+(0.)*np.sin(actin_angle_start)+(d_sp)*np.cos(actin_angle_start)
			elif pos_start-6<0.001:
				startx=actin_x_start+(15.)*np.cos(actin_angle_start)-(d_sp)*np.sin(actin_angle_start)
				starty=actin_y_start+(15.)*np.sin(actin_angle_start)+(d_sp)*np.cos(actin_angle_start)

		if actin_angle_end<=np.pi/4. or 7*np.pi/4.<actin_angle_end:
			if pos_end-1<0.001:
				endx=actin_x_end+(0.)*np.cos(actin_angle_end)-(d_sp)*np.sin(actin_angle_end)
				endy=actin_y_end+(0.)*np.sin(actin_angle_end)+(d_sp)*np.cos(actin_angle_end)
			elif pos_end-2<0.001:
				endx=actin_x_end+(-15.)*np.cos(actin_angle_end)-(d_sp)*np.sin(actin_angle_end)
				endy=actin_y_end+(-15.)*np.sin(actin_angle_end)+(d_sp)*np.cos(actin_angle_end)
			elif pos_end-3<0.001:
				endx=actin_x_end+(-15.)*np.cos(actin_angle_end)-(-d_sp)*np.sin(actin_angle_end)
				endy=actin_y_end+(-15.)*np.sin(actin_angle_end)+(-d_sp)*np.cos(actin_angle_end)
		elif (np.pi/4.<actin_angle_end and actin_angle_end<=3*np.pi/4.):
			if pos_end-1<0.001:
				endx=actin_x_end+(15.)*np.cos(actin_angle_end)-(d_sp)*np.sin(actin_angle_end)
				endy=actin_y_end+(15.)*np.sin(actin_angle_end)+(d_sp)*np.cos(actin_angle_end)
			elif pos_end-2<0.001:
				endx=actin_x_end+(0.)*np.cos(actin_angle_end)-(d_sp)*np.sin(actin_angle_end)
				endy=actin_y_end+(0.)*np.sin(actin_angle_end)+(d_sp)*np.cos(actin_angle_end)
			elif pos_end-3<0.001:
				endx=actin_x_end+(-15.)*np.cos(actin_angle_end)-(d_sp)*np.sin(actin_angle_end)
				endy=actin_y_end+(-15.)*np.sin(actin_angle_end)+(d_sp)*np.cos(actin_angle_end)
			#special case of anchor
			elif pos_end-5<0.001:
				endx=actin_x_end+(0.)*np.cos(actin_angle_end)-(-d_sp)*np.sin(actin_angle_end)
				endy=actin_y_end+(0.)*np.sin(actin_angle_end)+(-d_sp)*np.cos(actin_angle_end)
		elif (3*np.pi/4.<actin_angle_end and actin_angle_end<=5*np.pi/4.): 
			if pos_end-1<0.001:
				endx=actin_x_end+(15.)*np.cos(actin_angle_end)-(-d_sp)*np.sin(actin_angle_end)
				endy=actin_y_end+(15.)*np.sin(actin_angle_end)+(-d_sp)*np.cos(actin_angle_end)
			elif pos_end-2<0.001:
				endx=actin_x_end+(15.)*np.cos(actin_angle_end)-(d_sp)*np.sin(actin_angle_end)
				endy=actin_y_end+(15.)*np.sin(actin_angle_end)+(d_sp)*np.cos(actin_angle_end)
			elif pos_end-3<0.001:
				endx=actin_x_end+(0.)*np.cos(actin_angle_end)-(d_sp)*np.sin(actin_angle_end)
				endy=actin_y_end+(0.)*np.sin(actin_angle_end)+(d_sp)*np.cos(actin_angle_end)
		elif (5*np.pi/4.<actin_angle_end and actin_angle_end<=7*np.pi/4.): 
			if pos_end-1<0.001:
				endx=actin_x_end+(-15.)*np.cos(actin_angle_end)-(-d_sp)*np.sin(actin_angle_end)
				endy=actin_y_end+(-15.)*np.sin(actin_angle_end)+(-d_sp)*np.cos(actin_angle_end)
			elif pos_end-2<0.001:
				endx=actin_x_end+(0.)*np.cos(actin_angle_end)-(-d_sp)*np.sin(actin_angle_end)
				endy=actin_y_end+(0.)*np.sin(actin_angle_end)+(-d_sp)*np.cos(actin_angle_end)
			elif pos_end-3<0.001:
				endx=actin_x_end+(15.)*np.cos(actin_angle_end)-(-d_sp)*np.sin(actin_angle_end)
				endy=actin_y_end+(15.)*np.sin(actin_angle_end)+(-d_sp)*np.cos(actin_angle_end)

		init_spectrin_pos=np.zeros((39,3))

		d=np.sqrt(np.power(endy-starty,2)+np.power(endx-startx,2))
		alpha=np.arctan((endy-starty)/(endx-startx))
		delta_d=d/38.
		delta_d_small=(d-16*delta_d)/32.

		if d<20:
			z=z_plane+10
		else:
			z=z_plane				
				
		if np.random.rand()<0.5:
			for i in range(39):
				if i<14:
					init_spectrin_pos[i,0]=periodic_x(startx+i*delta_d*np.cos(alpha))
					init_spectrin_pos[i,1]=periodic_y(starty+i*delta_d*np.sin(alpha))
					init_spectrin_pos[i,2]=z+15*np.power(np.sin((i/15.*2)*np.pi*0.5),1)
				else:
					init_spectrin_pos[i,0]=periodic_x(startx+(15*delta_d+(i-15)*delta_d)*np.cos(alpha))
					init_spectrin_pos[i,1]=periodic_y(starty+(15*delta_d+(i-15)*delta_d)*np.sin(alpha))
					init_spectrin_pos[i,2]=z+15*np.sin(((i-12)/26.*2)*np.pi*0.5)
			
			spectrin = simulation.add_topology("spectrin", ["spectrinedge", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrinbinding", "spectrinbinding", "spectrinbinding", "spectrinbinding", "spectrinmidpoint", "spectrinmidpoint", "spectrincore", "spectrinbinding", "spectrincore", "spectrincore", "spectrincore", "spectrinbinding","spectrincore", "spectrincore", "spectrincore", "spectrinbinding", "spectrinbinding", "spectrinbinding", "spectrinbinding", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrinedge"], init_spectrin_pos)

			choices.append(1)
		else:
			for i in range(39):
				if i<=24:
					init_spectrin_pos[i,0]=periodic_x(startx+i*delta_d*np.cos(alpha))
					init_spectrin_pos[i,1]=periodic_y(starty+i*delta_d*np.sin(alpha))
					init_spectrin_pos[i,2]=z+15*np.sin(i/27.*2*np.pi*0.5)
				else:
					init_spectrin_pos[i,0]=periodic_x(startx+i*delta_d*np.cos(alpha))
					init_spectrin_pos[i,1]=periodic_y(starty+i*delta_d*np.sin(alpha))
					init_spectrin_pos[i,2]=z+15*np.sin((i-24)/14*2*np.pi*0.5)

			spectrin = simulation.add_topology("spectrin", ["spectrinedge", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrinbinding", "spectrinbinding", "spectrinbinding", "spectrinbinding", "spectrincore", "spectrincore", "spectrincore", "spectrinbinding", "spectrincore", "spectrincore","spectrincore", "spectrinbinding", "spectrincore","spectrinmidpoint", "spectrinmidpoint", "spectrinbinding", "spectrinbinding", "spectrinbinding", "spectrinbinding", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrinedge"], init_spectrin_pos)

			choices.append(2)

		for i in range(38):
			spectrin.get_graph().add_edge(i, i+1)


	#Initialize spectrin filaments    actin_x_start,actin_y_start,actin_x_end,actin_y_end,pos_start,pos_end,z,actin_angle_start,actin_angle_end

	choices=[]
	anchor_angle=np.pi/2.0

	#horizontals

	spectrin_initialization_periodic(x_coords[1],y_coords[1],x_coords[3],y_coords[1],5,2,z_pos,actin_angles[0],actin_angles[1],choices)
	spectrin_initialization_periodic(x_coords[3],y_coords[1],x_coords[5],y_coords[1],5,2,z_pos,actin_angles[1],actin_angles[2],choices)
	spectrin_initialization_periodic(x_coords[5],y_coords[1],x_coords[7],y_coords[1],5,2,z_pos,actin_angles[2],actin_angles[3],choices)
	spectrin_initialization_periodic(x_coords[7],y_coords[1],x_coords[9],y_coords[1],5,2,z_pos,actin_angles[3],actin_angles[4],choices)
	spectrin_initialization_periodic(x_coords[9],y_coords[1],x_coords[11],y_coords[1],5,2,z_pos,actin_angles[4],actin_angles[5],choices)
	spectrin_initialization_periodic(x_coords[11],y_coords[1],x_coords[13],y_coords[1],5,2,z_pos,actin_angles[5],actin_angles[6],choices)

	spectrin_initialization_periodic(x_coords[0],y_coords[2],x_coords[2],y_coords[2],5,2,z_pos,anchor_angle,actin_angles[7],choices)
	spectrin_initialization_periodic(x_coords[2],y_coords[2],x_coords[4],y_coords[2],5,2,z_pos,actin_angles[7],actin_angles[8],choices)
	spectrin_initialization_periodic(x_coords[4],y_coords[2],x_coords[6],y_coords[2],5,2,z_pos,actin_angles[8],actin_angles[9],choices)
	spectrin_initialization_periodic(x_coords[6],y_coords[2],x_coords[8],y_coords[2],5,2,z_pos,actin_angles[9],actin_angles[10],choices)
	spectrin_initialization_periodic(x_coords[8],y_coords[2],x_coords[10],y_coords[2],5,2,z_pos,actin_angles[10],actin_angles[11],choices)
	spectrin_initialization_periodic(x_coords[10],y_coords[2],x_coords[12],y_coords[2],5,2,z_pos,actin_angles[11],actin_angles[12],choices)
	spectrin_initialization_periodic(x_coords[12],y_coords[2],x_coords[14],y_coords[2],5,2,z_pos,actin_angles[12],anchor_angle,choices)

	spectrin_initialization_periodic(x_coords[1],y_coords[3],x_coords[3],y_coords[3],5,2,z_pos,actin_angles[13],actin_angles[14],choices)
	spectrin_initialization_periodic(x_coords[3],y_coords[3],x_coords[5],y_coords[3],5,2,z_pos,actin_angles[14],actin_angles[15],choices)
	spectrin_initialization_periodic(x_coords[5],y_coords[3],x_coords[7],y_coords[3],5,2,z_pos,actin_angles[15],actin_angles[16],choices)
	spectrin_initialization_periodic(x_coords[7],y_coords[3],x_coords[9],y_coords[3],5,2,z_pos,actin_angles[16],actin_angles[17],choices)
	spectrin_initialization_periodic(x_coords[9],y_coords[3],x_coords[11],y_coords[3],5,2,z_pos,actin_angles[17],actin_angles[18],choices)
	spectrin_initialization_periodic(x_coords[11],y_coords[3],x_coords[13],y_coords[3],5,2,z_pos,actin_angles[18],actin_angles[19],choices)

	spectrin_initialization_periodic(x_coords[0],y_coords[4],x_coords[2],y_coords[4],5,2,z_pos,anchor_angle,actin_angles[20],choices)
	spectrin_initialization_periodic(x_coords[2],y_coords[4],x_coords[4],y_coords[4],5,2,z_pos,actin_angles[20],actin_angles[21],choices)
	spectrin_initialization_periodic(x_coords[4],y_coords[4],x_coords[6],y_coords[4],5,2,z_pos,actin_angles[21],actin_angles[22],choices)
	spectrin_initialization_periodic(x_coords[6],y_coords[4],x_coords[8],y_coords[4],5,2,z_pos,actin_angles[22],actin_angles[23],choices)
	spectrin_initialization_periodic(x_coords[8],y_coords[4],x_coords[10],y_coords[4],5,2,z_pos,actin_angles[23],actin_angles[24],choices)
	spectrin_initialization_periodic(x_coords[10],y_coords[4],x_coords[12],y_coords[4],5,2,z_pos,actin_angles[24],actin_angles[25],choices)
	spectrin_initialization_periodic(x_coords[12],y_coords[4],x_coords[14],y_coords[4],5,2,z_pos,actin_angles[25],anchor_angle,choices)

	spectrin_initialization_periodic(x_coords[1],y_coords[5],x_coords[3],y_coords[5],5,2,z_pos,actin_angles[26],actin_angles[27],choices)
	spectrin_initialization_periodic(x_coords[3],y_coords[5],x_coords[5],y_coords[5],5,2,z_pos,actin_angles[27],actin_angles[28],choices)
	spectrin_initialization_periodic(x_coords[5],y_coords[5],x_coords[7],y_coords[5],5,2,z_pos,actin_angles[28],actin_angles[29],choices)
	spectrin_initialization_periodic(x_coords[7],y_coords[5],x_coords[9],y_coords[5],5,2,z_pos,actin_angles[29],actin_angles[30],choices)
	spectrin_initialization_periodic(x_coords[9],y_coords[5],x_coords[11],y_coords[5],5,2,z_pos,actin_angles[30],actin_angles[31],choices)
	spectrin_initialization_periodic(x_coords[11],y_coords[5],x_coords[13],y_coords[5],5,2,z_pos,actin_angles[31],actin_angles[32],choices)

	spectrin_initialization_periodic(x_coords[0],y_coords[6],x_coords[2],y_coords[6],5,2,z_pos,anchor_angle,actin_angles[33],choices)
	spectrin_initialization_periodic(x_coords[2],y_coords[6],x_coords[4],y_coords[6],5,2,z_pos,actin_angles[33],actin_angles[34],choices)
	spectrin_initialization_periodic(x_coords[4],y_coords[6],x_coords[6],y_coords[6],5,2,z_pos,actin_angles[34],actin_angles[35],choices)
	spectrin_initialization_periodic(x_coords[6],y_coords[6],x_coords[8],y_coords[6],5,2,z_pos,actin_angles[35],actin_angles[36],choices)
	spectrin_initialization_periodic(x_coords[8],y_coords[6],x_coords[10],y_coords[6],5,2,z_pos,actin_angles[36],actin_angles[37],choices)
	spectrin_initialization_periodic(x_coords[10],y_coords[6],x_coords[12],y_coords[6],5,2,z_pos,actin_angles[37],actin_angles[38],choices)
	spectrin_initialization_periodic(x_coords[12],y_coords[6],x_coords[14],y_coords[6],5,2,z_pos,actin_angles[38],anchor_angle,choices)

	spectrin_initialization_periodic(x_coords[1],y_coords[7],x_coords[3],y_coords[7],5,2,z_pos,actin_angles[39],actin_angles[40],choices)
	spectrin_initialization_periodic(x_coords[3],y_coords[7],x_coords[5],y_coords[7],5,2,z_pos,actin_angles[40],actin_angles[41],choices)
	spectrin_initialization_periodic(x_coords[5],y_coords[7],x_coords[7],y_coords[7],5,2,z_pos,actin_angles[41],actin_angles[42],choices)
	spectrin_initialization_periodic(x_coords[7],y_coords[7],x_coords[9],y_coords[7],5,2,z_pos,actin_angles[42],actin_angles[43],choices)
	spectrin_initialization_periodic(x_coords[9],y_coords[7],x_coords[11],y_coords[7],5,2,z_pos,actin_angles[43],actin_angles[44],choices)
	spectrin_initialization_periodic(x_coords[11],y_coords[7],x_coords[13],y_coords[7],5,2,z_pos,actin_angles[44],actin_angles[45],choices)

	#diagonals

	spectrin_initialization_periodic(x_coords[0],y_coords[0],x_coords[1],y_coords[1],5,3,z_pos,anchor_angle,actin_angles[0],choices)
	spectrin_initialization_periodic(x_coords[1],y_coords[1],x_coords[2],y_coords[0],6,2,z_pos,actin_angles[0],anchor_angle,choices)
	spectrin_initialization_periodic(x_coords[2],y_coords[0],x_coords[3],y_coords[1],5,3,z_pos,anchor_angle,actin_angles[1],choices)
	spectrin_initialization_periodic(x_coords[3],y_coords[1],x_coords[4],y_coords[0],6,2,z_pos,actin_angles[1],anchor_angle,choices)
	spectrin_initialization_periodic(x_coords[4],y_coords[0],x_coords[5],y_coords[1],5,3,z_pos,anchor_angle,actin_angles[2],choices)
	spectrin_initialization_periodic(x_coords[5],y_coords[1],x_coords[6],y_coords[0],6,2,z_pos,actin_angles[2],anchor_angle,choices)
	spectrin_initialization_periodic(x_coords[6],y_coords[0],x_coords[7],y_coords[1],5,3,z_pos,anchor_angle,actin_angles[3],choices)
	spectrin_initialization_periodic(x_coords[7],y_coords[1],x_coords[8],y_coords[0],6,2,z_pos,actin_angles[3],anchor_angle,choices)
	spectrin_initialization_periodic(x_coords[8],y_coords[0],x_coords[9],y_coords[1],5,3,z_pos,anchor_angle,actin_angles[4],choices)
	spectrin_initialization_periodic(x_coords[9],y_coords[1],x_coords[10],y_coords[0],6,2,z_pos,actin_angles[4],anchor_angle,choices)
	spectrin_initialization_periodic(x_coords[10],y_coords[0],x_coords[11],y_coords[1],5,3,z_pos,anchor_angle,actin_angles[5],choices)
	spectrin_initialization_periodic(x_coords[11],y_coords[1],x_coords[12],y_coords[0],6,2,z_pos,actin_angles[5],anchor_angle,choices)
	spectrin_initialization_periodic(x_coords[12],y_coords[0],x_coords[13],y_coords[1],5,3,z_pos,anchor_angle,actin_angles[6],choices)
	spectrin_initialization_periodic(x_coords[13],y_coords[1],x_coords[14],y_coords[0],6,2,z_pos,actin_angles[6],anchor_angle,choices)

	spectrin_initialization_periodic(x_coords[0],y_coords[2],x_coords[1],y_coords[1],5,1,z_pos,anchor_angle,actin_angles[0],choices)
	spectrin_initialization_periodic(x_coords[1],y_coords[1],x_coords[2],y_coords[2],4,3,z_pos,actin_angles[0],actin_angles[7],choices)
	spectrin_initialization_periodic(x_coords[2],y_coords[2],x_coords[3],y_coords[1],6,1,z_pos,actin_angles[7],actin_angles[1],choices)
	spectrin_initialization_periodic(x_coords[3],y_coords[1],x_coords[4],y_coords[2],4,3,z_pos,actin_angles[1],actin_angles[8],choices)
	spectrin_initialization_periodic(x_coords[4],y_coords[2],x_coords[5],y_coords[1],6,1,z_pos,actin_angles[8],actin_angles[2],choices)
	spectrin_initialization_periodic(x_coords[5],y_coords[1],x_coords[6],y_coords[2],4,3,z_pos,actin_angles[2],actin_angles[9],choices)
	spectrin_initialization_periodic(x_coords[6],y_coords[2],x_coords[7],y_coords[1],6,1,z_pos,actin_angles[9],actin_angles[3],choices)
	spectrin_initialization_periodic(x_coords[7],y_coords[1],x_coords[8],y_coords[2],4,3,z_pos,actin_angles[3],actin_angles[10],choices)
	spectrin_initialization_periodic(x_coords[8],y_coords[2],x_coords[9],y_coords[1],6,1,z_pos,actin_angles[10],actin_angles[4],choices)
	spectrin_initialization_periodic(x_coords[9],y_coords[1],x_coords[10],y_coords[2],4,3,z_pos,actin_angles[4],actin_angles[11],choices)
	spectrin_initialization_periodic(x_coords[10],y_coords[2],x_coords[11],y_coords[1],6,1,z_pos,actin_angles[11],actin_angles[5],choices)
	spectrin_initialization_periodic(x_coords[11],y_coords[1],x_coords[12],y_coords[2],4,3,z_pos,actin_angles[5],actin_angles[12],choices)
	spectrin_initialization_periodic(x_coords[12],y_coords[2],x_coords[13],y_coords[1],6,1,z_pos,actin_angles[12],actin_angles[6],choices)
	spectrin_initialization_periodic(x_coords[13],y_coords[1],x_coords[14],y_coords[2],4,2,z_pos,actin_angles[6],anchor_angle,choices)

	spectrin_initialization_periodic(x_coords[0],y_coords[2],x_coords[1],y_coords[3],5,3,z_pos,anchor_angle,actin_angles[13],choices)
	spectrin_initialization_periodic(x_coords[1],y_coords[3],x_coords[2],y_coords[2],6,1,z_pos,actin_angles[13],actin_angles[7],choices)
	spectrin_initialization_periodic(x_coords[2],y_coords[2],x_coords[3],y_coords[3],4,3,z_pos,actin_angles[7],actin_angles[14],choices)
	spectrin_initialization_periodic(x_coords[3],y_coords[3],x_coords[4],y_coords[2],6,1,z_pos,actin_angles[14],actin_angles[8],choices)
	spectrin_initialization_periodic(x_coords[4],y_coords[2],x_coords[5],y_coords[3],4,3,z_pos,actin_angles[8],actin_angles[15],choices)
	spectrin_initialization_periodic(x_coords[5],y_coords[3],x_coords[6],y_coords[2],6,1,z_pos,actin_angles[15],actin_angles[9],choices)
	spectrin_initialization_periodic(x_coords[6],y_coords[2],x_coords[7],y_coords[3],4,3,z_pos,actin_angles[9],actin_angles[16],choices)
	spectrin_initialization_periodic(x_coords[7],y_coords[3],x_coords[8],y_coords[2],6,1,z_pos,actin_angles[16],actin_angles[10],choices)
	spectrin_initialization_periodic(x_coords[8],y_coords[2],x_coords[9],y_coords[3],4,3,z_pos,actin_angles[10],actin_angles[17],choices)
	spectrin_initialization_periodic(x_coords[9],y_coords[3],x_coords[10],y_coords[2],6,1,z_pos,actin_angles[17],actin_angles[11],choices)
	spectrin_initialization_periodic(x_coords[10],y_coords[2],x_coords[11],y_coords[3],4,3,z_pos,actin_angles[11],actin_angles[18],choices)
	spectrin_initialization_periodic(x_coords[11],y_coords[3],x_coords[12],y_coords[2],6,1,z_pos,actin_angles[18],actin_angles[12],choices)
	spectrin_initialization_periodic(x_coords[12],y_coords[2],x_coords[13],y_coords[3],4,3,z_pos,actin_angles[12],actin_angles[19],choices)
	spectrin_initialization_periodic(x_coords[13],y_coords[3],x_coords[14],y_coords[2],6,2,z_pos,actin_angles[19],anchor_angle,choices)

	spectrin_initialization_periodic(x_coords[0],y_coords[4],x_coords[1],y_coords[3],5,1,z_pos,anchor_angle,actin_angles[13],choices)
	spectrin_initialization_periodic(x_coords[1],y_coords[3],x_coords[2],y_coords[4],4,3,z_pos,actin_angles[13],actin_angles[20],choices)
	spectrin_initialization_periodic(x_coords[2],y_coords[4],x_coords[3],y_coords[3],6,1,z_pos,actin_angles[20],actin_angles[14],choices)
	spectrin_initialization_periodic(x_coords[3],y_coords[3],x_coords[4],y_coords[4],4,3,z_pos,actin_angles[14],actin_angles[21],choices)
	spectrin_initialization_periodic(x_coords[4],y_coords[4],x_coords[5],y_coords[3],6,1,z_pos,actin_angles[21],actin_angles[15],choices)
	spectrin_initialization_periodic(x_coords[5],y_coords[3],x_coords[6],y_coords[4],4,3,z_pos,actin_angles[15],actin_angles[22],choices)
	spectrin_initialization_periodic(x_coords[6],y_coords[4],x_coords[7],y_coords[3],6,1,z_pos,actin_angles[22],actin_angles[16],choices)
	spectrin_initialization_periodic(x_coords[7],y_coords[3],x_coords[8],y_coords[4],4,3,z_pos,actin_angles[16],actin_angles[23],choices)
	spectrin_initialization_periodic(x_coords[8],y_coords[4],x_coords[9],y_coords[3],6,1,z_pos,actin_angles[23],actin_angles[17],choices)
	spectrin_initialization_periodic(x_coords[9],y_coords[3],x_coords[10],y_coords[4],4,3,z_pos,actin_angles[17],actin_angles[24],choices)
	spectrin_initialization_periodic(x_coords[10],y_coords[4],x_coords[11],y_coords[3],6,1,z_pos,actin_angles[24],actin_angles[18],choices)
	spectrin_initialization_periodic(x_coords[11],y_coords[3],x_coords[12],y_coords[4],4,3,z_pos,actin_angles[18],actin_angles[25],choices)
	spectrin_initialization_periodic(x_coords[12],y_coords[4],x_coords[13],y_coords[3],6,1,z_pos,actin_angles[25],actin_angles[19],choices)
	spectrin_initialization_periodic(x_coords[13],y_coords[3],x_coords[14],y_coords[4],4,2,z_pos,actin_angles[19],anchor_angle,choices)

	spectrin_initialization_periodic(x_coords[0],y_coords[4],x_coords[1],y_coords[5],5,3,z_pos,anchor_angle,actin_angles[26],choices)
	spectrin_initialization_periodic(x_coords[1],y_coords[5],x_coords[2],y_coords[4],6,1,z_pos,actin_angles[26],actin_angles[20],choices)
	spectrin_initialization_periodic(x_coords[2],y_coords[4],x_coords[3],y_coords[5],4,3,z_pos,actin_angles[20],actin_angles[27],choices)
	spectrin_initialization_periodic(x_coords[3],y_coords[5],x_coords[4],y_coords[4],6,1,z_pos,actin_angles[27],actin_angles[21],choices)
	spectrin_initialization_periodic(x_coords[4],y_coords[4],x_coords[5],y_coords[5],4,3,z_pos,actin_angles[21],actin_angles[28],choices)
	spectrin_initialization_periodic(x_coords[5],y_coords[5],x_coords[6],y_coords[4],6,1,z_pos,actin_angles[28],actin_angles[22],choices)
	spectrin_initialization_periodic(x_coords[6],y_coords[4],x_coords[7],y_coords[5],4,3,z_pos,actin_angles[22],actin_angles[29],choices)
	spectrin_initialization_periodic(x_coords[7],y_coords[5],x_coords[8],y_coords[4],6,1,z_pos,actin_angles[29],actin_angles[23],choices)
	spectrin_initialization_periodic(x_coords[8],y_coords[4],x_coords[9],y_coords[5],4,3,z_pos,actin_angles[23],actin_angles[30],choices)
	spectrin_initialization_periodic(x_coords[9],y_coords[5],x_coords[10],y_coords[4],6,1,z_pos,actin_angles[30],actin_angles[24],choices)
	spectrin_initialization_periodic(x_coords[10],y_coords[4],x_coords[11],y_coords[5],4,3,z_pos,actin_angles[24],actin_angles[31],choices)
	spectrin_initialization_periodic(x_coords[11],y_coords[5],x_coords[12],y_coords[4],6,1,z_pos,actin_angles[31],actin_angles[25],choices)
	spectrin_initialization_periodic(x_coords[12],y_coords[4],x_coords[13],y_coords[5],4,3,z_pos,actin_angles[25],actin_angles[32],choices)
	spectrin_initialization_periodic(x_coords[13],y_coords[5],x_coords[14],y_coords[4],6,2,z_pos,actin_angles[32],anchor_angle,choices)

	spectrin_initialization_periodic(x_coords[0],y_coords[6],x_coords[1],y_coords[5],5,1,z_pos,anchor_angle,actin_angles[26],choices)
	spectrin_initialization_periodic(x_coords[1],y_coords[5],x_coords[2],y_coords[6],4,3,z_pos,actin_angles[26],actin_angles[33],choices)
	spectrin_initialization_periodic(x_coords[2],y_coords[6],x_coords[3],y_coords[5],6,1,z_pos,actin_angles[33],actin_angles[27],choices)
	spectrin_initialization_periodic(x_coords[3],y_coords[5],x_coords[4],y_coords[6],4,3,z_pos,actin_angles[27],actin_angles[34],choices)
	spectrin_initialization_periodic(x_coords[4],y_coords[6],x_coords[5],y_coords[5],6,1,z_pos,actin_angles[34],actin_angles[28],choices)
	spectrin_initialization_periodic(x_coords[5],y_coords[5],x_coords[6],y_coords[6],4,3,z_pos,actin_angles[28],actin_angles[35],choices)
	spectrin_initialization_periodic(x_coords[6],y_coords[6],x_coords[7],y_coords[5],6,1,z_pos,actin_angles[35],actin_angles[29],choices)
	spectrin_initialization_periodic(x_coords[7],y_coords[5],x_coords[8],y_coords[6],4,3,z_pos,actin_angles[29],actin_angles[36],choices)
	spectrin_initialization_periodic(x_coords[8],y_coords[6],x_coords[9],y_coords[5],6,1,z_pos,actin_angles[36],actin_angles[30],choices)
	spectrin_initialization_periodic(x_coords[9],y_coords[5],x_coords[10],y_coords[6],4,3,z_pos,actin_angles[30],actin_angles[37],choices)
	spectrin_initialization_periodic(x_coords[10],y_coords[6],x_coords[11],y_coords[5],6,1,z_pos,actin_angles[37],actin_angles[31],choices)
	spectrin_initialization_periodic(x_coords[11],y_coords[5],x_coords[12],y_coords[6],4,3,z_pos,actin_angles[31],actin_angles[38],choices)
	spectrin_initialization_periodic(x_coords[12],y_coords[6],x_coords[13],y_coords[5],6,1,z_pos,actin_angles[38],actin_angles[32],choices)
	spectrin_initialization_periodic(x_coords[13],y_coords[5],x_coords[14],y_coords[6],4,2,z_pos,actin_angles[32],anchor_angle,choices)

	spectrin_initialization_periodic(x_coords[0],y_coords[6],x_coords[1],y_coords[7],5,3,z_pos,anchor_angle,actin_angles[39],choices)
	spectrin_initialization_periodic(x_coords[1],y_coords[7],x_coords[2],y_coords[6],6,1,z_pos,actin_angles[39],actin_angles[33],choices)
	spectrin_initialization_periodic(x_coords[2],y_coords[6],x_coords[3],y_coords[7],4,3,z_pos,actin_angles[33],actin_angles[40],choices)
	spectrin_initialization_periodic(x_coords[3],y_coords[7],x_coords[4],y_coords[6],6,1,z_pos,actin_angles[40],actin_angles[34],choices)
	spectrin_initialization_periodic(x_coords[4],y_coords[6],x_coords[5],y_coords[7],4,3,z_pos,actin_angles[34],actin_angles[41],choices)
	spectrin_initialization_periodic(x_coords[5],y_coords[7],x_coords[6],y_coords[6],6,1,z_pos,actin_angles[41],actin_angles[35],choices)
	spectrin_initialization_periodic(x_coords[6],y_coords[6],x_coords[7],y_coords[7],4,3,z_pos,actin_angles[35],actin_angles[42],choices)
	spectrin_initialization_periodic(x_coords[7],y_coords[7],x_coords[8],y_coords[6],6,1,z_pos,actin_angles[42],actin_angles[36],choices)
	spectrin_initialization_periodic(x_coords[8],y_coords[6],x_coords[9],y_coords[7],4,3,z_pos,actin_angles[36],actin_angles[43],choices)
	spectrin_initialization_periodic(x_coords[9],y_coords[7],x_coords[10],y_coords[6],6,1,z_pos,actin_angles[43],actin_angles[37],choices)
	spectrin_initialization_periodic(x_coords[10],y_coords[6],x_coords[11],y_coords[7],4,3,z_pos,actin_angles[37],actin_angles[44],choices)
	spectrin_initialization_periodic(x_coords[11],y_coords[7],x_coords[12],y_coords[6],6,1,z_pos,actin_angles[44],actin_angles[38],choices)
	spectrin_initialization_periodic(x_coords[12],y_coords[6],x_coords[13],y_coords[7],4,3,z_pos,actin_angles[38],actin_angles[45],choices)
	spectrin_initialization_periodic(x_coords[13],y_coords[7],x_coords[14],y_coords[6],6,2,z_pos,actin_angles[45],anchor_angle,choices)

	spectrin_initialization_periodic(x_coords[0],y_coords[8],x_coords[1],y_coords[7],5,1,z_pos,anchor_angle,actin_angles[39],choices)
	spectrin_initialization_periodic(x_coords[1],y_coords[7],x_coords[2],y_coords[8],4,2,z_pos,actin_angles[39],anchor_angle,choices)
	spectrin_initialization_periodic(x_coords[2],y_coords[8],x_coords[3],y_coords[7],5,1,z_pos,anchor_angle,actin_angles[40],choices)
	spectrin_initialization_periodic(x_coords[3],y_coords[7],x_coords[4],y_coords[8],4,2,z_pos,actin_angles[40],anchor_angle,choices)
	spectrin_initialization_periodic(x_coords[4],y_coords[8],x_coords[5],y_coords[7],5,1,z_pos,anchor_angle,actin_angles[41],choices)
	spectrin_initialization_periodic(x_coords[5],y_coords[7],x_coords[6],y_coords[8],4,2,z_pos,actin_angles[41],anchor_angle,choices)
	spectrin_initialization_periodic(x_coords[6],y_coords[8],x_coords[7],y_coords[7],5,1,z_pos,anchor_angle,actin_angles[42],choices)
	spectrin_initialization_periodic(x_coords[7],y_coords[7],x_coords[8],y_coords[8],4,2,z_pos,actin_angles[42],anchor_angle,choices)
	spectrin_initialization_periodic(x_coords[8],y_coords[8],x_coords[9],y_coords[7],5,1,z_pos,anchor_angle,actin_angles[43],choices)
	spectrin_initialization_periodic(x_coords[9],y_coords[7],x_coords[10],y_coords[8],4,2,z_pos,actin_angles[43],anchor_angle,choices)
	spectrin_initialization_periodic(x_coords[10],y_coords[8],x_coords[11],y_coords[7],5,1,z_pos,anchor_angle,actin_angles[44],choices)
	spectrin_initialization_periodic(x_coords[11],y_coords[7],x_coords[12],y_coords[8],4,2,z_pos,actin_angles[44],anchor_angle,choices)
	spectrin_initialization_periodic(x_coords[12],y_coords[8],x_coords[13],y_coords[7],5,1,z_pos,anchor_angle,actin_angles[45],choices)
	spectrin_initialization_periodic(x_coords[13],y_coords[7],x_coords[14],y_coords[8],4,2,z_pos,actin_angles[45],anchor_angle,choices)


	#positions_Gactin=random_sample(g_actins, -50, 50)
	#simulation.add_particles("Gactin", positions_Gactin)

	#positions_adducin = random_sample(adducins, -50, 50)
	#simulation.add_particles("adducin", positions_adducin)

	#positions_tropomodulin = random_sample(tropomodulins, -50, 50)
	#simulation.add_particles("tropomodulin", positions_tropomodulin)


	#positions_kahrp = random_sample(kahrps, -10, 50)
	#simulation.add_particles("kahrp", positions_kahrp)

	positions_anchor=np.zeros((22,3))

	positions_anchor[0][0]=x_coords[0]
	positions_anchor[0][1]=y_coords[0]
	positions_anchor[0][2]=z_pos
	positions_anchor[1][0]=x_coords[2]
	positions_anchor[1][1]=y_coords[0]
	positions_anchor[1][2]=z_pos
	positions_anchor[2][0]=x_coords[4]
	positions_anchor[2][1]=y_coords[0]
	positions_anchor[2][2]=z_pos
	positions_anchor[3][0]=x_coords[6]
	positions_anchor[3][1]=y_coords[0]
	positions_anchor[3][2]=z_pos
	positions_anchor[4][0]=x_coords[8]
	positions_anchor[4][1]=y_coords[0]
	positions_anchor[4][2]=z_pos
	positions_anchor[5][0]=x_coords[10]
	positions_anchor[5][1]=y_coords[0]
	positions_anchor[5][2]=z_pos
	positions_anchor[6][0]=x_coords[12]
	positions_anchor[6][1]=y_coords[0]
	positions_anchor[6][2]=z_pos
	positions_anchor[7][0]=x_coords[14]
	positions_anchor[7][1]=y_coords[0]
	positions_anchor[7][2]=z_pos

	positions_anchor[8][0]=x_coords[0]
	positions_anchor[8][1]=y_coords[2]
	positions_anchor[8][2]=z_pos
	positions_anchor[9][0]=x_coords[14]
	positions_anchor[9][1]=y_coords[2]
	positions_anchor[9][2]=z_pos

	positions_anchor[10][0]=x_coords[0]
	positions_anchor[10][1]=y_coords[4]
	positions_anchor[10][2]=z_pos
	positions_anchor[11][0]=x_coords[14]
	positions_anchor[11][1]=y_coords[4]
	positions_anchor[11][2]=z_pos

	positions_anchor[12][0]=x_coords[0]
	positions_anchor[12][1]=y_coords[6]
	positions_anchor[12][2]=z_pos
	positions_anchor[13][0]=x_coords[14]
	positions_anchor[13][1]=y_coords[6]
	positions_anchor[13][2]=z_pos

	positions_anchor[14][0]=x_coords[0]
	positions_anchor[14][1]=y_coords[8]
	positions_anchor[14][2]=z_pos
	positions_anchor[15][0]=x_coords[2]
	positions_anchor[15][1]=y_coords[8]
	positions_anchor[15][2]=z_pos
	positions_anchor[16][0]=x_coords[4]
	positions_anchor[16][1]=y_coords[8]
	positions_anchor[16][2]=z_pos
	positions_anchor[17][0]=x_coords[6]
	positions_anchor[17][1]=y_coords[8]
	positions_anchor[17][2]=z_pos
	positions_anchor[18][0]=x_coords[8]
	positions_anchor[18][1]=y_coords[8]
	positions_anchor[18][2]=z_pos
	positions_anchor[19][0]=x_coords[10]
	positions_anchor[19][1]=y_coords[8]
	positions_anchor[19][2]=z_pos
	positions_anchor[20][0]=x_coords[12]
	positions_anchor[20][1]=y_coords[8]
	positions_anchor[20][2]=z_pos
	positions_anchor[21][0]=x_coords[14]
	positions_anchor[21][1]=y_coords[8]
	positions_anchor[21][2]=z_pos

	#simulation.add_particles("anchor", positions_anchor)

	# In[12]:


	simulation.output_file = out_file
	if os.path.exists(simulation.output_file):
	    os.remove(simulation.output_file)
	simulation.observe.topologies(writeout)
	simulation.record_trajectory(writeout)
	simulation.observe.virial(stride=writeout,callback=lambda x: print(x))
	simulation.progress_output_stride = writeout

	#simulation.make_checkpoints(checkpoints, output_directory=cps, max_n_saves=2) #00

	# In[13]:

	def particles_callback(particles):
	    types, ids, positions = particles
	simulation.observe.particles(stride=writeout)
	simulation.observe.forces(stride=writeout,types=None)
	simulation.run(runtime, .01)

#run a follow up run (assumes a previous simulation to read the positions from)

def set_up_and_run_follow_up(x_base,y_base,out_file,cps,writeout,checkpoints,runtime,last_types,last_positions,system2,gamma):

	system2.periodic_boundary_conditions = [False, False, False]
	volume=0.001*x_base*0.001*y_base*0.1 #um^3  (constricted by potential in z-direction

	# In[4]:

	epsilon_k_k=system2.kbt*0.1

	epsilon_k_sp=epsilon_k_k
	epsilon_k_a=epsilon_k_k


	#definition of particles with diffusion constants

	system2.add_species("Gactin", 71.5 * readdy.units.um**2 / readdy.units.s)
	system2.add_species("adducin", 51.07 * readdy.units.um**2 / readdy.units.s)
	system2.add_species("tropomodulin", 29.59 * readdy.units.um**2 / readdy.units.s)
	system2.add_species("kahrp", 76.6 * readdy.units.um**2 / readdy.units.s) 

	system2.topologies.add_type("filament") 
	system2.add_topology_species("core", 1 * readdy.units.um**2 / readdy.units.s)  
	system2.add_topology_species("pointed", 1 * readdy.units.um**2 / readdy.units.s)   
	system2.add_topology_species("barbed", 0.53 * readdy.units.um**2 / readdy.units.s)  
	system2.add_topology_species("pointedcap", 1 * readdy.units.um**2 / readdy.units.s) 
	system2.add_topology_species("barbedcap", 0.53 * readdy.units.um**2 / readdy.units.s) 

	system2.topologies.add_type("spectrin")
	system2.add_topology_species("spectrincore", 41.1 * readdy.units.um**2 / readdy.units.s)
	system2.add_topology_species("spectrinedge", 41.1 * readdy.units.um**2 / readdy.units.s)
	system2.add_topology_species("spectrinmidpoint", 0.53 * readdy.units.um**2 / readdy.units.s)
	system2.add_topology_species("spectrinbinding", 41.1 * readdy.units.um**2 / readdy.units.s)

	system2.add_species("anchor", 0)

	# In[5]:


	#add box and plane potential

	for i in ["Gactin","adducin","kahrp","tropomodulin","core","pointed","barbed","pointedcap","barbedcap","spectrincore","spectrinbinding","spectrinedge","spectrinmidpoint","anchor"]: 
	    system2.potentials.add_box(
		particle_type=i, force_constant=100, origin=[-x_base*1.1, -y_base*1.1*0.5, -50], 
		extent=[2*x_base*1.1,y_base*1.1,100]
	    )
	for i in ["core","pointed","barbed","pointedcap","barbedcap","spectrinmidpoint"]:   
	    system2.potentials.add_box(
		particle_type=i, force_constant=10, origin=[-x_base*1.1, -y_base*1.1*0.5, -42], #confinement in plane in middle of box
		extent=[2*x_base*1.1,y_base*1.1, 10]
	    )


	# In[6]:


	particle_radii = {"Gactin": 3., "adducin": 4.2, "kahrp": 2.8, "tropomodulin": 7.25, "core": 3., "pointed": 3.,"pointedcap": 7.25, "barbed": 3., "barbedcap": 4.2, "spectrincore": 5.26/2.0, "spectrinedge": 5.26/2.0,"spectrinmidpoint":5.26/2.0,"spectrinbinding":5.26/2.0,"anchor": 5.0} #radius in nm

	epsilon_const=14*4.184

	sigma_sp_a=(particle_radii["spectrincore"]+particle_radii["Gactin"])#np.power(2,5/6.)*
	sigma_sp=(particle_radii["spectrincore"]+particle_radii["spectrincore"])
	sigma_a=(particle_radii["Gactin"]+particle_radii["Gactin"])
	sigma_k_a=(particle_radii["kahrp"]+particle_radii["Gactin"])
	sigma_k_sp=(particle_radii["kahrp"]+particle_radii["spectrincore"])
	sigma_k_k=(particle_radii["kahrp"]+particle_radii["kahrp"])

	bond_force_spectrin=36*np.power(2,2./3)*epsilon_const/np.power(sigma_sp,2)
	angle_force_spectrin=2.62
	repulsion_force_spectrin=36*np.power(2,2./3)*epsilon_const/np.power(sigma_sp,2)

	bond_force=36*np.power(2,2./3)*epsilon_const/np.power(sigma_a,2)
	angle_force=4280.
	repulsion_force=36*np.power(2,2./3)*epsilon_const/np.power(sigma_a,2)

	repulsion_force_sp_a=36*np.power(2,2./3)*epsilon_const/np.power(sigma_sp_a,2)
	repulsion_force_k_a=36*np.power(2,2./3)*epsilon_const/np.power(sigma_k_a,2)
	repulsion_force_k_sp=36*np.power(2,2./3)*epsilon_const/np.power(sigma_k_sp,2)



	#direction does not matter (AB=BA)
	system2.topologies.configure_harmonic_bond("core", "core", force_constant=bond_force, length=particle_radii["core"]+particle_radii["core"])
	system2.topologies.configure_harmonic_bond("core", "pointed", force_constant=bond_force, length=particle_radii["core"]+particle_radii["pointed"])
	system2.topologies.configure_harmonic_bond("core", "barbed", force_constant=bond_force, length=particle_radii["core"]+particle_radii["barbed"])
	system2.topologies.configure_harmonic_bond("core", "pointedcap", force_constant=bond_force, length=particle_radii["core"]+particle_radii["pointedcap"])
	system2.topologies.configure_harmonic_bond("core", "barbedcap", force_constant=bond_force, length=particle_radii["core"]+particle_radii["barbedcap"])
	system2.topologies.configure_harmonic_bond("pointed", "barbed", force_constant=bond_force, length=particle_radii["pointed"]+particle_radii["barbed"])
	system2.topologies.configure_harmonic_bond("pointedcap", "barbed", force_constant=bond_force, length=particle_radii["pointedcap"]+particle_radii["barbed"])
	system2.topologies.configure_harmonic_bond("pointed", "barbedcap", force_constant=bond_force, length=particle_radii["pointed"]+particle_radii["barbedcap"])

	system2.topologies.configure_harmonic_angle("core", "core", "core", force_constant=angle_force, equilibrium_angle=np.pi)
	system2.topologies.configure_harmonic_angle("pointed", "core", "core", force_constant=angle_force, equilibrium_angle=np.pi)
	system2.topologies.configure_harmonic_angle("barbed", "core", "core", force_constant=angle_force, equilibrium_angle=np.pi)
	system2.topologies.configure_harmonic_angle("pointedcap", "core", "core", force_constant=angle_force, equilibrium_angle=np.pi)
	system2.topologies.configure_harmonic_angle("barbedcap", "core", "core", force_constant=angle_force, equilibrium_angle=np.pi)
	system2.topologies.configure_harmonic_angle("pointed", "core", "barbed", force_constant=angle_force, equilibrium_angle=np.pi)
	system2.topologies.configure_harmonic_angle("barbed", "core", "pointed", force_constant=angle_force, equilibrium_angle=np.pi)
	system2.topologies.configure_harmonic_angle("pointedcap", "core", "barbed", force_constant=angle_force, equilibrium_angle=np.pi)
	system2.topologies.configure_harmonic_angle("barbedcap", "core", "pointed", force_constant=angle_force, equilibrium_angle=np.pi)

	system2.topologies.configure_harmonic_bond("spectrincore", "spectrincore", force_constant=bond_force_spectrin, length=particle_radii["spectrincore"]+particle_radii["spectrincore"])
	system2.topologies.configure_harmonic_bond("spectrincore", "spectrinedge", force_constant=bond_force_spectrin, length=particle_radii["spectrincore"]+particle_radii["spectrinedge"])
	system2.topologies.configure_harmonic_bond("spectrincore", "spectrinbinding", force_constant=bond_force_spectrin, length=particle_radii["spectrincore"]+particle_radii["spectrinbinding"])
	system2.topologies.configure_harmonic_bond("spectrincore", "spectrinmidpoint", force_constant=bond_force_spectrin, length=particle_radii["spectrincore"]+particle_radii["spectrinmidpoint"])
	system2.topologies.configure_harmonic_bond("spectrinbinding", "spectrinbinding", force_constant=bond_force_spectrin, length=particle_radii["spectrinbinding"]+particle_radii["spectrinbinding"])
	system2.topologies.configure_harmonic_bond("spectrinmidpoint", "spectrinbinding", force_constant=bond_force_spectrin, length=particle_radii["spectrinmidpoint"]+particle_radii["spectrinbinding"])
	system2.topologies.configure_harmonic_bond("spectrinmidpoint", "spectrinmidpoint", force_constant=bond_force_spectrin, length=particle_radii["spectrinmidpoint"]+particle_radii["spectrinmidpoint"])

	system2.topologies.configure_harmonic_angle("spectrincore", "spectrincore", "spectrincore", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system2.topologies.configure_harmonic_angle("spectrinedge", "spectrincore", "spectrincore", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system2.topologies.configure_harmonic_angle("spectrincore", "spectrinbinding", "spectrincore", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system2.topologies.configure_harmonic_angle("spectrincore", "spectrincore", "spectrinbinding", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system2.topologies.configure_harmonic_angle("spectrincore", "spectrinbinding", "spectrinbinding", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system2.topologies.configure_harmonic_angle("spectrinbinding", "spectrinbinding", "spectrinbinding", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system2.topologies.configure_harmonic_angle("spectrincore", "spectrincore", "spectrinmidpoint", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system2.topologies.configure_harmonic_angle("spectrincore", "spectrinmidpoint", "spectrincore", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system2.topologies.configure_harmonic_angle("spectrincore", "spectrinmidpoint", "spectrinmidpoint", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system2.topologies.configure_harmonic_angle("spectrinbinding", "spectrinmidpoint", "spectrinmidpoint", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system2.topologies.configure_harmonic_angle("spectrinbinding", "spectrinbinding", "spectrinmidpoint", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system2.topologies.configure_harmonic_angle("spectrinbinding", "spectrincore", "spectrinbinding", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system2.topologies.configure_harmonic_angle("spectrinbinding", "spectrincore", "spectrinmidpoint", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)
	system2.topologies.configure_harmonic_angle("spectrincore", "spectrinbinding", "spectrincore", force_constant=angle_force_spectrin, equilibrium_angle=np.pi)

	#harmonic repulsion of freely diffusing particles
	all_pairs = [("Gactin","Gactin"), ("adducin","adducin"), ("tropomodulin", "tropomodulin"), ("Gactin", "adducin"), ("Gactin", "tropomodulin"), ("adducin", "tropomodulin")]
	for pair in all_pairs:
	    distance = particle_radii[pair[0]] + particle_radii[pair[1]]
	    system2.potentials.add_harmonic_repulsion(pair[0], pair[1], force_constant=repulsion_force, 
		                                     interaction_distance=distance)

	#KAHRP with actin
	all_pairs = [("kahrp","Gactin"), ("kahrp", "adducin"), ("kahrp", "tropomodulin"), ("pointed", "kahrp"), ("pointedcap", "kahrp"), ("barbed", "kahrp"), ("barbedcap", "kahrp")]
	for pair in all_pairs:
	    distance = particle_radii[pair[0]] + particle_radii[pair[1]]
	    system2.potentials.add_harmonic_repulsion(pair[0], pair[1], force_constant=repulsion_force_k_a, 
		                                     interaction_distance=distance)

	#KAHRP with spectrin
	all_pairs = [("kahrp","spectrincore"),("kahrp","spectrinedge"), ("spectrinmidpoint","kahrp")] #spectrinbinding -> LJ
	for pair in all_pairs:
	    distance = particle_radii[pair[0]] + particle_radii[pair[1]]
	    system2.potentials.add_harmonic_repulsion(pair[0], pair[1], force_constant=repulsion_force_k_a, 
		                                     interaction_distance=distance)
	    
	#harmonic repulsion of freely diffusing particles with core
	all_pairs = [("core","core"),("Gactin","core"),("Gactin","pointedcap"),("Gactin","barbedcap"), ("adducin","core"), ("adducin","pointed"), ("adducin","pointedcap"), ("tropomodulin", "core"), ("tropomodulin", "barbed"), ("tropomodulin", "barbedcap"), ("tropomodulin", "pointedcap"), ("adducin", "barbedcap")]
	for pair in all_pairs:
	    distance = particle_radii[pair[0]] + particle_radii[pair[1]]
	    system2.potentials.add_harmonic_repulsion(pair[0], pair[1], force_constant=repulsion_force, 
		                                     interaction_distance=distance)

	#harmonic repulsion of actin filaments
	all_pairs = [("core","core"),("pointed","core"),("pointedcap","core"),("barbedcap","core"),("barbed","core"),("barbed","pointed"),("barbedcap","pointed"),("pointedcap","pointed"),("pointed","pointed"),("barbed","pointedcap"),("barbed","barbed"),("barbed","barbedcap"),("pointedcap","barbedcap"),("pointedcap","pointedcap"),("barbedcap","barbedcap")]
	for pair in all_pairs:
	    distance = particle_radii[pair[0]] + particle_radii[pair[1]]
	    system2.potentials.add_harmonic_repulsion(pair[0], pair[1], force_constant=repulsion_force, 
		                                     interaction_distance=distance)

	#harmonic repulsion of freely diffusing particles with spectrin
	all_pairs = [("Gactin","spectrincore"), ("adducin","spectrincore"), ("tropomodulin", "spectrincore"),("Gactin","spectrinbinding"), ("adducin","spectrinbinding"), ("tropomodulin", "spectrinbinding"),("Gactin","spectrinedge"), ("adducin","spectrinedge"), ("tropomodulin", "spectrinedge"),("Gactin","spectrinmidpoint"), ("adducin","spectrinmidpoint"), ("tropomodulin", "spectrinmidpoint")]
	for pair in all_pairs:
	    distance = particle_radii[pair[0]] + particle_radii[pair[1]]
	    system2.potentials.add_harmonic_repulsion(pair[0], pair[1], force_constant=repulsion_force_sp_a, 
		                                     interaction_distance=distance)
	    
	#harmonic repulsion of spectrin core 
	all_pairs = [("core","spectrincore"), ("pointed", "spectrincore"), ("pointedcap", "spectrincore"), ("barbed", "spectrincore"), ("barbedcap", "spectrincore")]
	for pair in all_pairs:
	    distance = particle_radii[pair[0]] + particle_radii[pair[1]]
	    system2.potentials.add_harmonic_repulsion(pair[0], pair[1], force_constant=repulsion_force_sp_a, 
		                                     interaction_distance=distance)

	#harmonic repulsion of spectrin heads
	all_pairs = [("spectrinedge","spectrinedge")]
	for pair in all_pairs:
	    distance = (particle_radii[pair[0]] + particle_radii[pair[1]])
	    system2.potentials.add_harmonic_repulsion(pair[0], pair[1], force_constant=repulsion_force_spectrin*2, 
		                                     interaction_distance=distance*2)

	#harmonic repulsion of spectrin midpoints
	all_pairs = [("spectrinmidpoint","core"),("spectrinmidpoint","barbed"),("spectrinmidpoint","barbedcap"),("spectrinmidpoint","pointed"),("spectrinmidpoint","pointedcap")]
	for pair in all_pairs:
	    distance = (particle_radii[pair[0]] + particle_radii[pair[1]])
	    system2.potentials.add_harmonic_repulsion(pair[0], pair[1], force_constant=repulsion_force_sp_a, 
		                                     interaction_distance=distance)

	all_pairs = [("spectrinbinding","core"),("spectrinbinding","barbed"),("spectrinbinding","barbedcap"), ("spectrinedge","barbed"), ("spectrinedge","barbedcap"),("spectrinbinding","pointed"),("spectrinedge","pointed"),("spectrinedge","pointedcap"),("spectrinbinding","pointedcap")]
	for pair in all_pairs:
	    distance = (particle_radii[pair[0]] + particle_radii[pair[1]])
	    system2.potentials.add_harmonic_repulsion(pair[0], pair[1], force_constant=repulsion_force_sp_a, 
		                                     interaction_distance=distance)

	all_pairs = [("spectrinbinding","spectrinbinding"),("spectrinbinding","spectrincore"),("spectrincore","spectrincore"),("spectrincore","spectrinedge"), ("spectrinmidpoint","spectrinmidpoint"), ("spectrinmidpoint","spectrinedge"), ("spectrinbinding","spectrinedge"), ("spectrinmidpoint","spectrincore"),("spectrinmidpoint","spectrinbinding")]
	for pair in all_pairs:
	    distance = (particle_radii[pair[0]] + particle_radii[pair[1]])
	    system2.potentials.add_harmonic_repulsion(pair[0], pair[1], force_constant=repulsion_force_spectrin, 
		                                     interaction_distance=distance)


	#Interaction of Spectrinedge with actin (core, pointed, barbed)

	system2.potentials.add_lennard_jones("spectrinedge","core", m=12, n=6, cutoff=20, shift=True, epsilon=epsilon_const, sigma=sigma_sp_a)
	system2.potentials.add_lennard_jones("spectrinedge","anchor", m=12, n=6, cutoff=20, shift=True, epsilon=epsilon_const*5, sigma=np.power(2,5/6.)*(particle_radii["spectrinedge"]+particle_radii["anchor"])/2.)


	#Binding to filament

	binding_rate_pointed=0.0000686  #in 1/ns  #0.0000068 
	binding_rate_barbed=0.000597 #in 1/ns  #0.00016
	binding_rate_tropomodulin=24324  #in 1/ns
	binding_rate_adducin=0.00329 #in 1/ns
	binding_radius=6.  #to middle of G-actin

	system2.topologies.add_spatial_reaction(
	    "BindPointed: filament(pointed) + (Gactin) -> filament(core--pointed)", 
	    rate=binding_rate_pointed, radius=binding_radius
	)
	system2.topologies.add_spatial_reaction(
	    "BindBarbed: filament(barbed) + (Gactin) -> filament(core--barbed)", 
	    rate=binding_rate_barbed, radius=binding_radius
	)
	system2.topologies.add_spatial_reaction(
	    "BindPointedCap: filament(pointed) + (tropomodulin) -> filament(core--pointedcap)", 
	    rate=binding_rate_tropomodulin, radius=binding_radius
	)
	system2.topologies.add_spatial_reaction(
	    "BindBarbedCap: filament(barbed) + (adducin) -> filament(core--barbedcap)", 
	    rate=binding_rate_adducin, radius=binding_radius
	)


	# In[8]:

	system2.topologies.add_structural_reaction("pointed_dissociation","filament", pointed_dissociation, pointed_rate)
	system2.topologies.add_structural_reaction("barbed_dissociation", "filament", barbed_dissociation, barbed_rate)
	system2.topologies.add_structural_reaction("pointedcap_dissociation", "filament", pointedcap_dissociation, pointedcap_rate)
	system2.topologies.add_structural_reaction("barbedcap_dissociation", "filament", barbedcap_dissociation, barbedcap_rate)

	# In[9]:


	simulation2 = system2.simulation(kernel="CPU")


	# In[10]: Add particles

	initialize_particles(last_positions,last_types,simulation2,gamma)
	
	# In[12]:


	simulation2.output_file = out_file
	if os.path.exists(simulation2.output_file):
	    os.remove(simulation2.output_file)
	#simulation2.observe.topologies(writeout)
	simulation2.record_trajectory(writeout)
	simulation2.observe.virial(stride=writeout,callback=lambda x: print(x))
	simulation2.progress_output_stride = writeout

	#simulation2.make_checkpoints(checkpoints, output_directory=cps, max_n_saves=2) #00

	# In[13]:

	def particles_callback(particles):
	    types, ids, positions = particles
	    #print("Particle 5 has type {}, id {}, and position {}."
	    #      .format(types[5], ids[5], positions[5])
	simulation2.observe.particles(stride=writeout)
	simulation2.observe.forces(stride=writeout,types=None)

	simulation2.run(runtime, .01)

#function to apply a gradual strain to the network by shifting the beads accordingly

def transform(old_pos,gamma):
	new_pos=np.zeros(np.shape(old_pos))
	for i in range(len(old_pos)):
		p=old_pos[i]
		y_scaling=(p[1]+y_base*0.5)/y_base
		new_pos[i][0]=p[0]+y_scaling*gamma
		new_pos[i][1]=p[1]
		new_pos[i][2]=p[2]
	#print("old",old_pos)
	#print("new",new_pos)
	return new_pos

#functions to extract the final positions from the similation data and set up the cytoskeleton accordingly for the next simulation 

def initialize_actin_new(init_top_pos,pc,bc,sim):
	if pc:
        	struc=["pointedcap"]
	else:
        	struc=["pointed"]
        
	for i in range(len(init_top_pos)-2):
        	struc.append("core")
        
	if bc:
        	struc.append("barbedcap")
	else:
        	struc.append("barbed")
        
	top = sim.add_topology("filament", struc,init_top_pos)
	for i in range(len(init_top_pos)-1):
        	top.get_graph().add_edge(i, i+1)
    
    
def initialize_particles(positions,type_names,simu,gamma):        
	print("gamma",gamma)

	#46 actins
	ind=0
	for i in range(46):

		while type_names[ind]!='pointedcap' and type_names[ind]!='pointed':
			if type_names[ind]=='adducin':
				position_adducin=np.array([positions[ind]])
				simu.add_particles("adducin", position_adducin)
				#ind+=1
			elif type_names[ind]=='tropomodulin':
				position_tropomodulin=np.array([positions[ind]])
				simu.add_particles("tropomodulin", position_tropomodulin)			
				#ind+=1
			else:
				print("type error")
			ind+=1

		if type_names[ind]=='pointedcap':
            		pc=True
		elif type_names[ind]=='pointed':
            		pc=False
		s=ind
		ind+=1
		while type_names[ind]=='core':
            		ind+=1
		if type_names[ind]=='barbedcap':
            		bc=True
		elif type_names[ind]=='barbed':
            		bc=False
		e=ind+1
		initialize_actin_new(transform(positions[s:e],gamma),pc,bc,simu)
		ind+=1

	while type_names[ind]!='spectrinedge':
		if type_names[ind]=='adducin':
			position_adducin=np.array([positions[ind]])
			simu.add_particles("adducin", position_adducin)
			#ind+=1
		elif type_names[ind]=='tropomodulin':
			position_tropomodulin=np.array([positions[ind]])
			simu.add_particles("tropomodulin", position_tropomodulin)			
			#ind+=1
		else:
			print("type error")
		ind+=1
		

	#157 spectrins
	for i in range(157):
		struc_sp=[]
		s=ind
		e=ind+39
		init_spectrin_pos=transform(positions[s:e],gamma)
		for j in range(39):
			struc_sp.append(type_names[ind])
			ind+=1
		spectrin2=simu.add_topology("spectrin", struc_sp, init_spectrin_pos)

		for i in range(38):
			spectrin2.get_graph().add_edge(i, i+1)

	#free particles
	sG=ind
	e=len(type_names)

	if sG<e:
		while type_names[ind]=="Gactin":
			ind+=1
		eG=ind
		if sG<eG:
			positions_Gactin=positions[sG:eG]

			simu.add_particles("Gactin", positions_Gactin)

	sa=ind
	
	if sa<e:
		while type_names[ind]=="adducin":
			ind+=1
		ea=ind
		if sa<ea:
			positions_adducin=positions[sa:ea]
			simu.add_particles("adducin", positions_adducin)

	st=ind

	if st<e:
		while type_names[ind]=="tropomodulin":
			ind+=1
		et=ind
		if st<et:
			positions_tropomodulin=positions[st:et]
			simu.add_particles("tropomodulin", positions_tropomodulin)



###########################################################################################################################################	
#set up a complete shearing sequence

n_shear=200
gamma=np.linspace(0.005,1,n_shear)
rs=['_' + str(int(round(i*1000))) for i in gamma]

#lattice constant of hexagonal network (value changed from 60 to 120)
a=60
h=np.sqrt(0.75)*a
x_base=7*a
y_base=8*h

runtime_first=2000000
writeout_first=int(runtime_first/100)
checkpoints_first=int(runtime_first/10)

runtime=int(1600)
writeout=int(runtime/100)
checkpoints=int(runtime/10)

#first run

mean_xy=[]
mean_yx=[]

set_up_and_run_first(x_base,y_base,'stress_follow_up.h5',"checkpoints_stress_follow_up/",writeout_first,checkpoints_first,runtime_first)

traj = readdy.Trajectory('stress_follow_up.h5')
traj.convert_to_xyz(particle_radii = {"Gactin": 3., "adducin": 4.2, "kahrp": 2.8, "tropomodulin": 7.25, "core": 3., "pointed": 3.,"pointedcap": 7.25, "barbed": 3., "barbedcap": 4.2, "spectrincore": 5.26/2.0, "spectrinedge": 5.26/2.0,"spectrinmidpoint":5.26/2.0,"spectrinbinding":5.26/2.0,"anchor": 5.0})

extract_quantities_own(traj,mean_xy,mean_yx)

t, types, ids, pos = traj.read_observable_particles()
type_names=[]
for t_id in types[-1]:
    	type_names.append(traj.species_name(t_id))
positions=pos[-1]

x=x_base
y=y_base

#following runs

n_shear=200
gamma=np.linspace(0.005,1,n_shear)
rs=['_' + str(int(round(i*1000))) for i in gamma]

sf=[i*x_base for i in gamma]
sf_last=0.0

for i in range(n_shear):

	print(str(i) + " of 200 shear steps done.")
	
	r=rs[i]
	stretch_factor=sf[i]

	system_new = readdy.ReactionDiffusionSystem(box_size=[2*x*1.2, y*1.2, 100+100])
	set_up_and_run_follow_up(x,y,'stress_follow_up_stretch.h5','checkpoints_stress_follow_up_stretch/',writeout,checkpoints,runtime,type_names,positions,system_new,stretch_factor-sf_last)

	traj2 = readdy.Trajectory('stress_follow_up_stretch.h5')
	sf_last=stretch_factor

	extract_quantities_own(traj2,mean_xy,mean_yx)

	t, types, ids, pos = traj2.read_observable_particles()
	type_names=[]
	for t_id in types[-1]:
	    	type_names.append(traj.species_name(t_id))
	positions=pos[-1]

#save stress data to file
np.save("mean_xy.npy",mean_xy)
np.save("mean_yx.npy",mean_yx)



