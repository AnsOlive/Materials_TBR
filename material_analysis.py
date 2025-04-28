import openmc
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
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
initial_breeder_widths = np.linspace(0, 25, 6)

def run_openmc_sim(enrich, width):
    ##################################################################################################
    # Main Simulation
    ##################################################################################################

    ##################################################################################################
    # Materials
    ##################################################################################################
   
    # Defines the natural berillium reflector
    Breeder = openmc.Material(name="Breeding Layer")
    Breeder.set_density("g/cm3", 1.845)
    Breeder.add_element("Be", 100)
    materials = openmc.Materials([Breeder])

    # Creates all the salts you want to define below at a constant density of 1.91 g/cm^3, 
    # so they are the same as FLiBe, expcept for ones with LiPb which have a density closer to 8.11 g/cm^3
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
        salt.set_density('g/cm3', 8.11)
        salt.add_element('Li', 40, enrichment=enrich, enrichment_target='Li6')
        salt.add_element('Cu', 20)
        salt.add_element('Pb', 40)

    elif salt_name == 'LiGaPb':
        salt = openmc.Material(name="LiGaPb")
        salt.set_density('g/cm3', 8.11)
        salt.add_element('Li', 35, enrichment=enrich, enrichment_target='Li6')
        salt.add_element('Ga', 10)
        salt.add_element('Pb', 55)

    elif salt_name == 'LiSrPb':
        salt = openmc.Material(name="LiSrPb")
        salt.set_density('g/cm3', 8.11)
        salt.add_element('Li', 30, enrichment=enrich, enrichment_target='Li6')
        salt.add_element('Sr', 50)
        salt.add_element('Pb', 20)

    elif salt_name == 'LiPbZn':
        salt = openmc.Material(name="LiPbZn")
        salt.set_density('g/cm3', 8.11)
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
    reflector_radius = 600 + width
    boundary = outer_radius + 1

    # Define all boundaries
    sphere_inner = openmc.Sphere(r=inner_radius)
    sphere_reflector = openmc.Sphere(r=reflector_radius)
    sphere_outer = openmc.Sphere(r=outer_radius)
    sphere_leakage = openmc.Sphere(r=boundary, boundary_type='reflective')

    # Define all regions
    center_void = -sphere_inner
    region_reflector = +sphere_inner & -sphere_reflector
    region_salt = +sphere_reflector & -sphere_outer
    outer_void = +sphere_outer & -sphere_leakage

    # Define all cells
    center_cell = openmc.Cell(fill=None, region=center_void)
    cell_reflector = openmc.Cell(fill=Breeder, region=region_reflector)
    cell_salt = openmc.Cell(fill=salt, region=region_salt)
    outer_cell = openmc.Cell(fill=None, region=outer_void)

    # Create the universe
    universe = openmc.Universe(cells=[center_cell, cell_salt, cell_reflector, outer_cell])
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

# Run all the initial enrichments and widths
coarse_results = []
for width in initial_breeder_widths:
    for enrich in initial_enrichments:
        TBR_mean, TBR_std_dev = run_openmc_sim(enrich, width)
        coarse_results.append([width, enrich, TBR_mean, TBR_std_dev])

# Convert results into a NumPy array for structured processing
coarse_results = np.array(coarse_results)
coarse_widths = coarse_results[:, 0]
coarse_enrichments = coarse_results[:, 1]
coarse_TBR_values = coarse_results[:, 2]
max_index = np.argmax(coarse_TBR_values)
coarse_optimal_width = coarse_widths[max_index]
coarse_optimal_enrichment = coarse_enrichments[max_index]

##################################################################################################
# Iterative Refinement
##################################################################################################

# Refine the results around the optimal values
enrich_threshold = 1.0
width_threshold = 1.0

# Define the ranges for the refined search
near_maximum_enrichments = np.arange(max(10, coarse_optimal_enrichment - 5), min(100, coarse_optimal_enrichment + 5), enrich_threshold)
near_maximum_breeder_thicknesses = np.arange(max(0, coarse_optimal_width - 2.5), min(25, coarse_optimal_width + 2.5), width_threshold)

# Run the OpenMC simulation for each enrichment in the refined range
refined_results = []
for width in near_maximum_breeder_thicknesses:
    for enrich in near_maximum_enrichments:
        TBR_mean, TBR_std_dev = run_openmc_sim(enrich, width)
        refined_results.append([width, enrich, TBR_mean, TBR_std_dev])

# Combine the results from the initial run and the local refinement
refined_results = np.array(refined_results)

##################################################################################################
# Postprocessing
##################################################################################################

# Combine the results from the coarse and refined searches
combined_results = np.vstack([coarse_results, refined_results])

# Extract values for widths and enrichments
widths = combined_results[:, 0]
enrichments = combined_results[:, 1]
TBR_values = combined_results[:, 2]

# Create a grid for surface plotting
X_unique = np.unique(widths)
Y_unique = np.unique(enrichments)
X, Y = np.meshgrid(X_unique, Y_unique)

# Interpolate w values onto the grid
Z = griddata((widths, enrichments), TBR_values, (X, Y), method='linear')

# Find optimal values
max_index = np.nanargmax(TBR_values)
optimal_width = widths[max_index]
optimal_enrichment = enrichments[max_index]
max_tbr = TBR_values[max_index]

##################################################################################################
# Total Plots
##################################################################################################

# Heatmap Plot
fig1, ax1 = plt.subplots(figsize=(10, 7))
heatmap = ax1.contourf(X, Y, Z, levels=20, cmap="viridis")
cbar = plt.colorbar(heatmap)
ax1.set_xlabel("Breeder Blanket Width (cm)")
ax1.set_ylabel("Lithium Enrichment (%)")
ax1.set_title("TBR vs. Width vs. Enrichment (Heatmap)")
cbar.set_label("TBR Value")
plt.savefig(f"tbr_heatmap_{salt_name}.png", dpi=300)
plt.close(fig1)

# 3D Plot
fig2 = plt.figure(figsize=(10, 7))
ax2 = fig2.add_subplot(111, projection='3d')
ax2.plot_surface(X, Y, Z, cmap="viridis")
ax2.set_xlabel("Breeder Blanket Width (cm)")
ax2.set_ylabel("Lithium Enrichment (%)")
ax2.set_zlabel("TBR")
ax2.set_title("TBR vs. Width vs. Enrichment (3D Surface Plot)")
plt.savefig(f"tbr_3d_plot_{salt_name}.png", dpi=300)
plt.close(fig2)

##################################################################################################
# Coarse Plots
##################################################################################################

# TBR vs. Width at Optimal Enrichment for Coarse Search
mask_enrichment = np.isclose(coarse_enrichments, coarse_optimal_enrichment, atol=1.0)
widths_at_opt_enrich = coarse_widths[mask_enrichment]
tbr_at_opt_enrich = coarse_TBR_values[mask_enrichment]
sorted_indices = np.argsort(widths_at_opt_enrich)
widths_sorted = widths_at_opt_enrich[sorted_indices]
tbr_sorted = tbr_at_opt_enrich[sorted_indices]

fig3, ax3 = plt.subplots(figsize=(10, 6))
ax3.plot(widths_sorted, tbr_sorted, linestyle='-')
ax3.set_xlabel("Breeder Blanket Width (cm)")
ax3.set_ylabel("TBR")
ax3.set_title(f"TBR vs Width at {coarse_optimal_enrichment:.2f}% Enrichment (Initial Data) for {salt_name}")
ax3.set_xlim(0, 25)
plt.grid(True)
plt.savefig(f"tbr_vs_width_at_enrich_{salt_name}.png", dpi=300)
plt.close(fig3)

# TBR vs. Enrichment at Optimal Width for Coarse Search
mask_width = np.isclose(coarse_widths, coarse_optimal_width, atol=1.0)
enrichments_at_opt_width = coarse_enrichments[mask_width]
tbr_at_opt_width = coarse_TBR_values[mask_width]
sorted_indices = np.argsort(enrichments_at_opt_width)
enrichments_sorted = enrichments_at_opt_width[sorted_indices]
tbr_sorted_width = tbr_at_opt_width[sorted_indices]

fig4, ax4 = plt.subplots(figsize=(10, 6))
ax4.plot(enrichments_sorted, tbr_sorted_width, linestyle='--')
ax4.set_xlabel("Lithium Enrichment (%)")
ax4.set_ylabel("TBR")
ax4.set_title(f"TBR vs Enrichment at {coarse_optimal_width:.2f} cm Width (Initial Data) for {salt_name}")
ax4.grid(True)
ax4.set_xlim(10, 100)
plt.savefig(f"tbr_vs_enrichment_at_width_coarse_{salt_name}.png", dpi=300)
plt.close(fig4)

##################################################################################################
# Save Data
##################################################################################################

# Print results in a readable format
with open(f'raw_data_{salt_name}.txt', 'w') as file:
    file.write("### Coarse Search Data (initial) ###\n")
    file.write("Width (cm) | Enrichment (%) | TBR Mean | TBR Std Dev\n")
    for width, enrichment, tbr_mean, tbr_std in coarse_results:
        file.write(f"{width:10.2f} | {enrichment:14.2f} | {tbr_mean:8.5f} | {tbr_std:10.5f}\n")

    file.write("\n### Refined Search Data (local refinement) ###\n")
    file.write("Width (cm) | Enrichment (%) | TBR Mean | TBR Std Dev\n")
    for width, enrichment, tbr_mean, tbr_std in refined_results:
        file.write(f"{width:10.2f} | {enrichment:14.2f} | {tbr_mean:8.5f} | {tbr_std:10.5f}\n")

    file.write("\n### TBR Summary Statistics:\n")
    file.write(f"Maximum TBR for {salt_name}: {max_tbr:.5f}\n")
    file.write(f"Optimal breeder width: {optimal_width:.2f} cm\n")
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
