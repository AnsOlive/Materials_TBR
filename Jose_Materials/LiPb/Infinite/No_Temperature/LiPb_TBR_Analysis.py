# Joseangel Lima-Tapia
# Research Project Spring 2025
# Started: 3/18/2025

import openmc
from openmc.stats import Point, Discrete, Isotropic

####################
##### MATERIALS ####
####################

LiPb = openmc.Material(name = 'Lithium Lead')
LiPb.add_element('Pb' , 1)
LiPb.add_element('Li' , 1)
LiPb.set_density('g/cm3' , 10.218)
#LiPb.temperature = (873 + 923) / 2

Materials = openmc.Materials([LiPb])
Materials.export_to_xml()

####################
##### GEOMETRY #####
####################

hollow_sphere = openmc.Sphere(r = 10)
inf_sphere = openmc.Sphere(r = 20000 , boundary_type = 'reflective')

Void_Cell = openmc.Cell(region = -hollow_sphere)
LiPb_Cell = openmc.Cell(region = +hollow_sphere & -inf_sphere , fill = LiPb)

Geometry = openmc.Geometry([LiPb_Cell, Void_Cell])
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
Tally.filters = [openmc.CellFilter(LiPb_Cell)]

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
#Settings.temperature = {
#    'tolerance': 135,  # K
#    'method': 'interpolation'
#}

Settings.export_to_xml()

####################
#### Simulation ####
####################

openmc.run()