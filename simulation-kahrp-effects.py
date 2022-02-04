#script to set up and run a simulation of a periodic cytosceleton patch with KAHRP particles present
#used to generate the data displayed in figure 5 and 6
#author: Julia Jaeger

import numpy as np
import matplotlib.pyplot as plt
import readdy
import os

def simulation(factor_k_k,factor_k_spectrin,factor_k_actin,factor_k_ankyrin,out_file): #factor is epsilon in units of kbT, 'parameter-scan-01.h5'
	x_stretch=1.0
	y_stretch=1.0

	a=88
	h=np.sqrt(0.75)*a

	x_base=2*a
	y_base=4*h

	z_pos=-40

	system = readdy.ReactionDiffusionSystem(box_size=[x_base*x_stretch, y_base*y_stretch, 100+100])
	system.periodic_boundary_conditions = [True, True, False]
	volume=0.001*x_base*x_stretch*0.001*y_base*y_stretch*0.1 #um^3  (constricted by potential in z-direction


	g_actins = 20#126
	kahrps=200
	tropomodulins= 5
	adducins= 5


	# In[4]:

	epsilon_k_k=system.kbt*factor_k_k

	epsilon_k_spectrin=system.kbt*factor_k_spectrin
	epsilon_k_actin=system.kbt*factor_k_actin
	epsilon_k_ankyrin=system.kbt*factor_k_ankyrin


	#definition of particles with diffusion constants

	system.add_species("Gactin", 71.5 * readdy.units.um**2 / readdy.units.s)
	system.add_species("adducin", 51.07 * readdy.units.um**2 / readdy.units.s)
	system.add_species("tropomodulin", 29.59 * readdy.units.um**2 / readdy.units.s)
	system.add_species("kahrp", 76.6 * readdy.units.um**2 / readdy.units.s) 

	system.topologies.add_type("filament") 
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


	# In[5]:


	#add box and plane potential

	for i in ["Gactin","adducin","kahrp","tropomodulin","core","pointed","barbed","pointedcap","barbedcap","spectrincore","spectrinbinding","spectrinedge","spectrinmidpoint"]: 
	    system.potentials.add_box(
		particle_type=i, force_constant=100, origin=[-3000, -3000, -50], 
		extent=[6000,6000,100]
	    )
	for i in ["core","pointed","barbed","pointedcap","barbedcap","spectrinmidpoint"]:   
	    system.potentials.add_box(
		particle_type=i, force_constant=10, origin=[-3000, -3000, -42], #confinement in plane in middle of box
		extent=[6000, 6000, 4]
	    )


	# In[6]:


	particle_radii = {"Gactin": 3., "adducin": 4.2, "kahrp": 2.8, "tropomodulin": 7.25, "core": 3., "pointed": 3.,"pointedcap": 7.25, "barbed": 3., "barbedcap": 4.2, "spectrincore": 5.26/2.0, "spectrinedge": 5.26/2.0,"spectrinmidpoint":5.26/2.0,"spectrinbinding":5.26/2.0} #radius in nm

	epsilon_const=14*4.184

	sigma_sp_a=(particle_radii["spectrincore"]+particle_radii["Gactin"])
	sigma_sp=(particle_radii["spectrincore"]+particle_radii["spectrincore"])
	sigma_a=(particle_radii["Gactin"]+particle_radii["Gactin"])
	sigma_k_a=(particle_radii["kahrp"]+particle_radii["Gactin"])
	sigma_k_sp=(particle_radii["kahrp"]+particle_radii["spectrincore"])
	sigma_k_k=(particle_radii["kahrp"]+particle_radii["kahrp"])

	bond_force_spectrin=36*np.power(2,2./3)*epsilon_const/np.power(sigma_sp,2)
	angle_force_spectrin=4.28
	repulsion_force_spectrin=36*np.power(2,2./3)*epsilon_const/np.power(sigma_sp,2)

	bond_force=36*np.power(2,2./3)*epsilon_const/np.power(sigma_a,2)
	angle_force=4280.
	repulsion_force=36*np.power(2,2./3)*epsilon_const/np.power(sigma_a,2)

	repulsion_force_sp_a=36*np.power(2,2./3)*epsilon_const/np.power(sigma_sp_a,2)
	repulsion_force_k_a=36*np.power(2,2./3)*epsilon_const/np.power(sigma_k_a,2)
	repulsion_force_k_sp=36*np.power(2,2./3)*epsilon_const/np.power(sigma_k_sp,2)



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
	    system.potentials.add_harmonic_repulsion(pair[0], pair[1], force_constant=repulsion_force_k_a,interaction_distance=distance)

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



	#binding of spectrin to actin
	system.potentials.add_lennard_jones("spectrinedge","core", m=12, n=6, cutoff=20, shift=True, epsilon=epsilon_const, sigma=sigma_sp_a)

	#kahrp binding
	system.potentials.add_lennard_jones("spectrinbinding","kahrp", m=12, n=6, cutoff=20, shift=True, epsilon=epsilon_k_spectrin, sigma=sigma_k_sp)
	system.potentials.add_lennard_jones("kahrp","kahrp", m=12, n=6, cutoff=20, shift=True, epsilon=epsilon_k_k, sigma=sigma_k_k)
	system.potentials.add_lennard_jones("kahrp","core", m=12, n=6, cutoff=20, shift=True, epsilon=epsilon_k_actin, sigma=sigma_k_a)
	system.potentials.add_lennard_jones("kahrp","spectrinmidpoint", m=12, n=6, cutoff=20, shift=True, epsilon=epsilon_k_ankyrin, sigma=sigma_k_sp)


	# In[7]:

	#Binding to filament

	binding_rate_pointed=0.0000686  #in 1/ns   
	binding_rate_barbed=0.000597 #in 1/ns  
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


	# In[8]:


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

	system.topologies.add_structural_reaction("pointed_dissociation","filament", pointed_dissociation, pointed_rate)
	system.topologies.add_structural_reaction("barbed_dissociation", "filament", barbed_dissociation, barbed_rate)
	system.topologies.add_structural_reaction("pointedcap_dissociation", "filament", pointedcap_dissociation, pointedcap_rate)
	system.topologies.add_structural_reaction("barbedcap_dissociation", "filament", barbedcap_dissociation, barbedcap_rate)

	# In[9]:


	simulation = system.simulation(kernel="CPU")


	# In[10]:

	x_coords=[-x_base*0.5,-x_base*0.5+a*0.5,-x_base*0.5+a,-x_base*0.5+a*1.5,x_base*0.5]
	y_coords=[-y_base*0.5,-y_base*0.5+h,-y_base*0.5+2*h,-y_base*0.5+3*h,y_base*0.5]

	x_coords=[x_stretch*i for i in x_coords]
	y_coords=[y_stretch*i for i in y_coords]

	actin_angles=np.random.rand(8)*2*np.pi


	#initialize actin filaments
	def periodic_x(x):
		if x>=x_coords[0] and x<=x_coords[4]:
			x_new=x
		elif x>x_coords[4]:
			x_new=x_coords[0]+x-x_coords[4]
		elif x<x_coords[4]:
			x_new=x_coords[4]+(x-x_coords[0])
		return x_new

	def periodic_y(y):
		if y>=y_coords[0] and y<=y_coords[4]:
			y_new=y
		elif y>y_coords[4]:
			y_new=y_coords[0]+y-y_coords[4]
		elif y<y_coords[4]:
			y_new=y_coords[4]+(y-y_coords[0])
		return y_new

	def initialize_actin(x,y,z,alpha):
		init_top_pos = np.array([
		    [periodic_x(x+(-15.)*np.cos(alpha)),periodic_y(y+(-15.)*np.sin(alpha)),z],
		    [periodic_x(x+(-9.)*np.cos(alpha)),periodic_y(y+(-9.)*np.sin(alpha)),z],
		    [periodic_x(x+(-3.)*np.cos(alpha)),periodic_y(y+(-3.)*np.sin(alpha)),z],
		    [periodic_x(x+(3.)*np.cos(alpha)),periodic_y(y+(3.)*np.sin(alpha)),z],
		    [periodic_x(x+(9.)*np.cos(alpha)),periodic_y(y+(9.)*np.sin(alpha)),z],
		    [periodic_x(x+(15.)*np.cos(alpha)),periodic_y(y+(15.)*np.sin(alpha)),z]
		])
		top = simulation.add_topology("filament", ["pointedcap", "core", "core", "core", "core", "barbedcap"], init_top_pos)
		top.get_graph().add_edge(0, 1)
		top.get_graph().add_edge(1, 2)
		top.get_graph().add_edge(2, 3)
		top.get_graph().add_edge(3, 4)
		top.get_graph().add_edge(4, 5)	



	initialize_actin(x_coords[1],y_coords[0],z_pos,actin_angles[0])
	initialize_actin(x_coords[3],y_coords[0],z_pos,actin_angles[1])

	initialize_actin(x_coords[2],y_coords[1],z_pos,actin_angles[2])
	initialize_actin(x_coords[4],y_coords[1],z_pos,actin_angles[3])

	initialize_actin(x_coords[1],y_coords[2],z_pos,actin_angles[4])
	initialize_actin(x_coords[3],y_coords[2],z_pos,actin_angles[5])

	initialize_actin(x_coords[2],y_coords[3],z_pos,actin_angles[6])
	initialize_actin(x_coords[4],y_coords[3],z_pos,actin_angles[7])


	def spectrin_initialization_periodic(actin_x_start,actin_y_start,actin_x_end,actin_y_end,pos_start,pos_end,z,actin_angle_start,actin_angle_end):
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
			if pos_start-4<0.001:
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

		if abs(endy-starty)>2*h and endy<-h:
			endy=y_coords[4]+(endy-y_coords[0])
		if abs(endy-starty)>2*h and starty<-h:
			starty=y_coords[4]+(starty-y_coords[0])
		if endx<0 and startx>0:
			endx=x_coords[4]+(endx-x_coords[0])

		d=np.sqrt(np.power(endy-starty,2)+np.power(endx-startx,2))
		alpha=np.arctan((endy-starty)/(endx-startx))
		delta_d=d/38.
		delta_d_small=(d-16*delta_d)/32.


		
				
				
		if np.random.rand()<0.5:
			for i in range(39):
				if i<14:
					init_spectrin_pos[i,0]=periodic_x(startx+i*delta_d*np.cos(alpha))
					init_spectrin_pos[i,1]=periodic_y(starty+i*delta_d*np.sin(alpha))
					init_spectrin_pos[i,2]=z+25*np.sin(i/13.*2*np.pi*0.5)
				elif i<=21:
					init_spectrin_pos[i,0]=periodic_x(startx+(12*delta_d+(i-12)*delta_d_small*2.)*np.cos(alpha))
					init_spectrin_pos[i,1]=periodic_y(starty+(12*delta_d+(i-12)*delta_d_small*2.)*np.sin(alpha))
					init_spectrin_pos[i,2]=z-8*np.power(np.sin((i-14)/7.*2*np.pi*0.5),2)
				else:
					init_spectrin_pos[i,0]=periodic_x(startx+(12*delta_d+9*delta_d_small*2+(i-21)*delta_d_small)*np.cos(alpha))
					init_spectrin_pos[i,1]=periodic_y(starty+(12*delta_d+9*delta_d_small*2+(i-21)*delta_d_small)*np.sin(alpha))
					init_spectrin_pos[i,2]=z+25*np.sin((i-21)/17.*2*np.pi*0.5)
			
			spectrin = simulation.add_topology("spectrin", ["spectrinedge", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrinbinding", "spectrinbinding", "spectrinbinding", "spectrinbinding", "spectrinmidpoint", "spectrinmidpoint", "spectrincore", "spectrinbinding", "spectrincore", "spectrincore", "spectrincore", "spectrinbinding","spectrincore", "spectrincore", "spectrincore", "spectrinbinding", "spectrinbinding", "spectrinbinding", "spectrinbinding", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrinedge"], init_spectrin_pos)
		else:
			for i in range(39):
				if i<16:
					init_spectrin_pos[i,0]=periodic_x(startx+i*delta_d_small*np.cos(alpha))
					init_spectrin_pos[i,1]=periodic_y(starty+i*delta_d_small*np.sin(alpha))
					init_spectrin_pos[i,2]=z+25*np.sin(i/16.*2*np.pi*0.5)
				elif i<=24:
					init_spectrin_pos[i,0]=periodic_x(startx+(15*delta_d_small+(i-15)*delta_d_small*2.)*np.cos(alpha))
					init_spectrin_pos[i,1]=periodic_y(starty+(15*delta_d_small+(i-15)*delta_d_small*2.)*np.sin(alpha))
					init_spectrin_pos[i,2]=z-8*np.sin(i/8.*2*np.pi*0.5)
				else:
					init_spectrin_pos[i,0]=periodic_x(startx+i*delta_d*np.cos(alpha))
					init_spectrin_pos[i,1]=periodic_y(starty+i*delta_d*np.sin(alpha))
					init_spectrin_pos[i,2]=z+25*np.power(np.sin((i-25)/13.*2*np.pi*0.5),2)

			spectrin = simulation.add_topology("spectrin", ["spectrinedge", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrinbinding", "spectrinbinding", "spectrinbinding", "spectrinbinding", "spectrincore", "spectrincore", "spectrincore", "spectrinbinding", "spectrincore", "spectrincore","spectrincore", "spectrinbinding", "spectrincore","spectrinmidpoint", "spectrinmidpoint", "spectrinbinding", "spectrinbinding", "spectrinbinding", "spectrinbinding", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrincore", "spectrinedge"], init_spectrin_pos)

		for i in range(38):
			spectrin.get_graph().add_edge(i, i+1)


	#Initialize spectrin filaments    actin_x_start,actin_y_start,actin_x_end,actin_y_end,pos_start,pos_end,z,actin_angle_start,actin_angle_end

	spectrin_initialization_periodic(x_coords[1],y_coords[0],x_coords[3],y_coords[0],5,2,z_pos,actin_angles[0],actin_angles[1])
	spectrin_initialization_periodic(x_coords[3],y_coords[0],x_coords[1],y_coords[0],5,2,z_pos,actin_angles[1],actin_angles[0])

	spectrin_initialization_periodic(x_coords[0],y_coords[1],x_coords[2],y_coords[1],5,2,z_pos,actin_angles[3],actin_angles[2])
	spectrin_initialization_periodic(x_coords[2],y_coords[1],x_coords[4],y_coords[1],5,2,z_pos,actin_angles[2],actin_angles[3])

	spectrin_initialization_periodic(x_coords[1],y_coords[2],x_coords[3],y_coords[2],5,2,z_pos,actin_angles[4],actin_angles[5])
	spectrin_initialization_periodic(x_coords[3],y_coords[2],x_coords[1],y_coords[2],5,2,z_pos,actin_angles[5],actin_angles[4])

	spectrin_initialization_periodic(x_coords[0],y_coords[3],x_coords[2],y_coords[3],5,2,z_pos,actin_angles[7],actin_angles[6])
	spectrin_initialization_periodic(x_coords[2],y_coords[3],x_coords[4],y_coords[3],5,2,z_pos,actin_angles[6],actin_angles[7])

	spectrin_initialization_periodic(x_coords[0],y_coords[1],x_coords[1],y_coords[0],6,1,z_pos,actin_angles[3],actin_angles[0])
	spectrin_initialization_periodic(x_coords[1],y_coords[0],x_coords[2],y_coords[1],4,3,z_pos,actin_angles[0],actin_angles[2])
	spectrin_initialization_periodic(x_coords[2],y_coords[1],x_coords[3],y_coords[0],6,1,z_pos,actin_angles[2],actin_angles[1])
	spectrin_initialization_periodic(x_coords[3],y_coords[0],x_coords[4],y_coords[1],4,3,z_pos,actin_angles[1],actin_angles[3])

	spectrin_initialization_periodic(x_coords[0],y_coords[1],x_coords[1],y_coords[2],4,3,z_pos,actin_angles[3],actin_angles[4])
	spectrin_initialization_periodic(x_coords[1],y_coords[2],x_coords[2],y_coords[1],6,1,z_pos,actin_angles[4],actin_angles[2])
	spectrin_initialization_periodic(x_coords[2],y_coords[1],x_coords[3],y_coords[2],4,3,z_pos,actin_angles[2],actin_angles[5])
	spectrin_initialization_periodic(x_coords[3],y_coords[2],x_coords[4],y_coords[1],6,1,z_pos,actin_angles[5],actin_angles[3])

	spectrin_initialization_periodic(x_coords[0],y_coords[3],x_coords[1],y_coords[2],6,1,z_pos,actin_angles[7],actin_angles[4])
	spectrin_initialization_periodic(x_coords[1],y_coords[2],x_coords[2],y_coords[3],4,3,z_pos,actin_angles[4],actin_angles[6])
	spectrin_initialization_periodic(x_coords[2],y_coords[3],x_coords[3],y_coords[2],6,1,z_pos,actin_angles[6],actin_angles[5])
	spectrin_initialization_periodic(x_coords[3],y_coords[2],x_coords[4],y_coords[3],4,3,z_pos,actin_angles[5],actin_angles[7])

	spectrin_initialization_periodic(x_coords[0],y_coords[3],x_coords[1],y_coords[4],4,3,z_pos,actin_angles[7],actin_angles[0])
	spectrin_initialization_periodic(x_coords[1],y_coords[4],x_coords[2],y_coords[3],6,1,z_pos,actin_angles[0],actin_angles[6])
	spectrin_initialization_periodic(x_coords[2],y_coords[3],x_coords[3],y_coords[4],4,3,z_pos,actin_angles[6],actin_angles[1])
	spectrin_initialization_periodic(x_coords[3],y_coords[4],x_coords[4],y_coords[3],6,1,z_pos,actin_angles[1],actin_angles[7])


	# In[11]:

	def random_sample(N, z_min, z_max):
		positions=np.zeros((N,3))
		for i in range(N):
			positions[i][0]=np.random.rand()*x_base*x_stretch - x_base*x_stretch* 0.5
			positions[i][1]=np.random.rand()*y_base*y_stretch - y_base*y_stretch* 0.5
			positions[i][2]=z_min+np.random.rand()*(z_max-z_min)
		return positions


	positions_Gactin=random_sample(g_actins, -50, 50)
	simulation.add_particles("Gactin", positions_Gactin)

	positions_adducin = random_sample(adducins, -50, 50)
	simulation.add_particles("adducin", positions_adducin)

	positions_tropomodulin = random_sample(tropomodulins, -50, 50)
	simulation.add_particles("tropomodulin", positions_tropomodulin)


	positions_kahrp = random_sample(kahrps, -10, 50)
	simulation.add_particles("kahrp", positions_kahrp)



	# In[12]:


	simulation.output_file = out_file
	if os.path.exists(simulation.output_file):
		os.remove(simulation.output_file)
	simulation.observe.particle_positions(stride=10000,types=["kahrp","spectrinbinding","spectrinmidpoint","spectrinedge","core"])
	simulation.observe.particles(stride=10000)
	simulation.record_trajectory(10000)
	simulation.progress_output_stride=10000

	simulation.run(20000000, .01)  #10000000

	

