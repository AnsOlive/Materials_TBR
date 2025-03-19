# Joseangel Lima-Tapia
# Research Project Spring 2025
# Started: 3/5/2025

import openmc
from openmc.stats import Point, Discrete, Isotropic

####################
##### MATERIALS ####
####################

FLiBe = openmc.Material(name = 'FLiBe')
FLiBe.add_element('F' , .573) # From HYLIFE paper
FLiBe.add_element('Li' , .282)
FLiBe.add_element('Be' , .145)
FLiBe.add_s_alpha_beta('c_Be' , fraction = .145)
FLiBe.set_density('g/cm3' , 1.94)
FLiBe.temperature = (823 + 923) / 2

Materials = openmc.Materials([FLiBe])
Materials.export_to_xml()

####################
##### GEOMETRY #####
####################

hollow_sphere = openmc.Sphere(r = 10)
inf_sphere = openmc.Sphere(r = 20000 , boundary_type = 'reflective')

Void_Cell = openmc.Cell(region = -hollow_sphere)
FliBe_Cell = openmc.Cell(region = +hollow_sphere & -inf_sphere , fill = FLiBe)

Geometry = openmc.Geometry([FliBe_Cell, Void_Cell])
Geometry.export_to_xml()

####################
###### Source ######
####################

Source = openmc.IndependentSource(                 # OpenMC already assumes the source is at the origin but its good practice to define the position
    space = Point(xyz = (0 , 0 , 0)),              # Centered on the origin
    angle = Isotropic(),                           # Isotropic particle distrubution
    energy = Discrete([14.1e6] , [1.0]),           # Discrete is saying there's a 14.1 MeV energy distrubution
    particle = 'neutron'                           # I believe this is naturally set but better practice to define it
)

####################
###### Tallies #####
####################

Tally = openmc.Tally(name = 'Tritium Production')
Tally.scores = ['(n,Xt)']
Tally.filters = [openmc.CellFilter(FliBe_Cell)]

Tallies = openmc.Tallies()
Tallies.append(Tally)

Tallies.export_to_xml()

####################
##### Settings #####
####################

Settings = openmc.Settings()

Settings.source = Source           # Defines my central neutron for the simulation
Settings.batches = 100             # 100 total rounds
Settings.particles = 10000         # more particles = lower uncertainty
Settings.run_mode = 'fixed source' # changes from k (eigenvalue) to source simulations
Settings.temperature = {
    'tolerance': 50,  # K
    'method': 'interpolation'
}

Settings.export_to_xml()

####################
#### Simulation ####
####################

openmc.run()