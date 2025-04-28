import openmc
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

##################################################################################################
# Initial Setup
##################################################################################################

# Find the salt name from command line arguments
parser = argparse.ArgumentParser(description='Run OpenMC simulation with specified salt.')
parser.add_argument('--salt', type=str, required=True, help='Name of the salt material to use')
args = parser.parse_args()
salt_name = args.salt

# Set the initial enrichment values
initial_enrichments = np.linspace(10, 100, 10)

def run_openmc_sim(enrich):
    ##################################################################################################
    # Main Simulation
    ##################################################################################################

    ##################################################################################################
    # Materials
    ##################################################################################################
   
    # Initializes materials 
    materials = openmc.Materials()

    # Creates though all the salts you want to define below at a constant density of 1.91 g/cm^3, so they are the same as FLiBe
    # No density specified in the paper and theorical densities are difficult to know

    if salt_name == 'FLiBe':
        salt = openmc.Material(name='FLiBe')
        salt.set_density('g/cm3', 1.91)
        salt.add_element('F', 57.3)
        salt.add_element('Li', 28.2, enrichment=enrich, enrichment_target='Li6')
        salt.add_element('Be', 14.5)

    elif salt_name == 'LiBaBi':
        salt = openmc.Material(name="LiBaBi")
        salt.set_density('g/cm3', 1.91)
        salt.add_element('Li', 20, enrichment=enrich, enrichment_target='Li6')
        salt.add_element('Ba', 10)
        salt.add_element('Bi', 70)

    elif salt_name == 'LiPbBa':
        salt = openmc.Material(name="LiPbBa")
        salt.set_density('g/cm3', 1.91)
        salt.add_element('Li', 25, enrichment=enrich, enrichment_target='Li6')
        salt.add_element('Pb', 60)
        salt.add_element('Ba', 15)

    elif salt_name == 'LiSnZn':
        salt = openmc.Material(name="LiSnZn")
        salt.set_density('g/cm3', 1.91)
        salt.add_element('Li', 65, enrichment=enrich, enrichment_target='Li6')
        salt.add_element('Sn', 25)
        salt.add_element('Zn', 10)

    elif salt_name == 'LiCuPb':
        salt = openmc.Material(name="LiCuPb")
        salt.set_density('g/cm3', 1.91)
        salt.add_element('Li', 40, enrichment=enrich, enrichment_target='Li6')
        salt.add_element('Cu', 20)
        salt.add_element('Pb', 40)

    elif salt_name == 'LiGaPb':
        salt = openmc.Material(name="LiGaPb")
        salt.set_density('g/cm3', 1.91)
        salt.add_element('Li', 35, enrichment=enrich, enrichment_target='Li6')
        salt.add_element('Ga', 10)
        salt.add_element('Pb', 55)

    elif salt_name == 'LiSrPb':
        salt = openmc.Material(name="LiSrPb")
        salt.set_density('g/cm3', 1.91)
        salt.add_element('Li', 30, enrichment=enrich, enrichment_target='Li6')
        salt.add_element('Sr', 50)
        salt.add_element('Pb', 20)

    elif salt_name == 'LiPbZn':
        salt = openmc.Material(name="LiPbZn")
        salt.set_density('g/cm3', 1.91)
        salt.add_element('Li', 30, enrichment=enrich, enrichment_target='Li6')
        salt.add_element('Pb', 60)
        salt.add_element('Zn', 10)

    elif salt_name == 'LiNaSn':
        salt = openmc.Material(name="LiNaSn")
        salt.set_density('g/cm3', 1.91)
        salt.add_element('Li', 55, enrichment=enrich, enrichment_target='Li6')
        salt.add_element('Na', 30)
        salt.add_element('Sn', 15)

    # Adds the selected salt to the simulation
    materials.append(salt)
    materials.export_to_xml()

    ##################################################################################################
    # Geometry
    ##################################################################################################

    # Defines all the radii for this iteration
    inner_radius, outer_radius = 600, 200000
    boundary = outer_radius + 1

    # Define all boundaries
    sphere_inner = openmc.Sphere(r=inner_radius)
    sphere_outer = openmc.Sphere(r=outer_radius)
    sphere_leakage = openmc.Sphere(r=boundary, boundary_type='reflective')

    # Define all regions
    center_void = -sphere_inner
    region_salt = +sphere_inner & -sphere_outer
    outer_void = +sphere_outer & -sphere_leakage

    # Define all cells
    center_cell = openmc.Cell(fill=None, region=center_void)
    cell_salt = openmc.Cell(fill=salt, region=region_salt)
    outer_cell = openmc.Cell(fill=None, region=outer_void)

    # Create the universe
    universe = openmc.Universe(cells=[center_cell, cell_salt, outer_cell])
    geometry = openmc.Geometry(universe)

    # Export the geometry to XML
    geometry.export_to_xml()

    ##################################################################################################
    # Settings
    ##################################################################################################

    # Define the source
    source = openmc.IndependentSource()
    source.particle = 'neutron'
    source.energy = openmc.stats.Discrete([14.1e6], [1.0])
    source.angle = openmc.stats.Isotropic()
    source.space = openmc.stats.Point((0.0, 0.0, 0.0))

    # Simulation settings
    settings = openmc.Settings()
    settings.batches = 50
    settings.inactive = 0
    settings.particles = 500
    settings.run_mode = 'fixed source'
    settings.source = source
    settings.export_to_xml()
        
    ##################################################################################################
    # Tallies
    ##################################################################################################

    # Initialize the tallies
    tallies = openmc.Tallies()

    # Looks at the overall TBR
    TBR_tally = openmc.Tally(name='TBR')
    TBR_tally.filters = [openmc.MaterialFilter([salt])]
    TBR_tally.scores = ['(n,Xt)']
    tallies.append(TBR_tally)

    tallies.export_to_xml()

    ##################################################################################################
    # Run the OpenMC Simulation
    ##################################################################################################
    openmc.run()

    ##################################################################################################
    # Add Results
    ##################################################################################################
    # Use the written results to add to the results array
    with openmc.StatePoint(f"statepoint.{settings.batches}.h5") as statepoint:
        TBR_mean = statepoint.get_tally(name="TBR").mean.item()
        TBR_std_dev = statepoint.get_tally(name="TBR").std_dev.item()

    # Output the results
    return TBR_mean, TBR_std_dev

##################################################################################################
# Initial Run
##################################################################################################

# Run all the initial enrichments
coarse_results = []
for enrich in initial_enrichments:
    TBR_mean, TBR_std_dev = run_openmc_sim(enrich)
    coarse_results.append([enrich, TBR_mean, TBR_std_dev])

# Convert results into a NumPy array for structured processing
coarse_results = np.array(coarse_results)
max_index = np.argmax(coarse_results[:, 1])
optimal_enrichment = coarse_results[max_index, 0]
max_tbr = coarse_results[max_index, 1]

##################################################################################################
# Iterative Refinement
##################################################################################################

# Refine the enrichment range around the optimal enrichment
enrich_threshold = 1.0

# Define the range of enrichments to explore around the optimal enrichment
enrichments = np.arange(max(10, optimal_enrichment - 5), min(100, optimal_enrichment + 5) + 0.5, enrich_threshold)

# Run the OpenMC simulation for each enrichment in the refined range
fine_results = []
for enrich in enrichments:
    TBR_mean, TBR_std_dev = run_openmc_sim(enrich)
    fine_results.append([enrich, TBR_mean, TBR_std_dev])
fine_results = np.array(fine_results)

# Combine the results from the initial run and the local refinement
combined_results = np.vstack([coarse_results, fine_results])
combined_results = np.array(combined_results)

##################################################################################################
# Plots
##################################################################################################

# Extract width, enrichment, and TBR values which are needed to generate plots
enrichments = combined_results[:, 0]
TBR_values = combined_results[:, 1]

# Find maximum TBR and corresponding enrichment for all enrichments
max_index = np.argmax(combined_results[:, 1])
optimal_enrichment = combined_results[max_index, 0]
max_tbr = combined_results[max_index, 1]

# Sort by enrichment before plotting
sorted_indices = np.argsort(enrichments)
enrichments = enrichments[sorted_indices]
TBR_values = TBR_values[sorted_indices]

# Plot TBR vs Enrichment
fig4, ax4 = plt.subplots(figsize=(10, 6))
ax4.plot(enrichments, TBR_values, linestyle='--')
ax4.set_xlabel("Lithium Enrichment (%)")
ax4.set_ylabel("TBR")
ax4.set_title(f"TBR vs Enrichment for {salt_name}")
plt.grid(True)
plt.savefig(f"tbr_vs_enrichment_{salt_name}.png", dpi=300)
plt.close(fig4)

##################################################################################################
# Save Data
##################################################################################################

# Print results
with open(f'raw_data_{salt_name}.txt', 'w') as file:
    file.write("### Coarse Search Data (initial) ###\n")
    file.write("Enrichment (%) | TBR Mean | TBR Std Dev\n")
    for enrichment, tbr_mean, tbr_std in coarse_results:
        file.write(f"{enrichment:14.2f} | {tbr_mean:8.5f} | {tbr_std:10.5f}\n")

    file.write("\n### Refined Search Data (local refinement) ###\n")
    file.write("Enrichment (%) | TBR Mean | TBR Std Dev\n")
    for enrichment, tbr_mean, tbr_std in fine_results:
        file.write(f"{enrichment:14.2f} | {tbr_mean:8.5f} | {tbr_std:10.5f}\n")

    file.write("\n### TBR Summary Statistics:\n")
    file.write(f"Maximum TBR for {salt_name}: {max_tbr:.5f}\n")
    file.write(f"Optimal lithium enrichment: {optimal_enrichment:.2f} %\n")
    file.write(f"\nMean TBR: {np.mean(TBR_values):.5f}\n")
    file.write(f"Median TBR: {np.median(TBR_values):.5f}\n")

print(f"Maximum TBR for {salt_name}: {max_tbr:.5f}")
print(f"Optimal breeder width: {optimal_width:.2f} cm")
print(f"Optimal lithium enrichment: {optimal_enrichment:.2f} %")

# Clean up OpenMC temporary files
for file in ['geometry.xml', 'materials.xml', 'settings.xml', 'summary.h5', f'statepoint.{settings.batches}.h5', 'tallies.xml']:
    if os.path.exists(file):
        os.remove(file)
