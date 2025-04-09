import openmc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import argparse
import os
import glob

##################################################################################################
# Initial Setup
##################################################################################################

# Parse arguments to get needed data from batch script
parser = argparse.ArgumentParser(description='Run OpenMC simulation with specified salt.')
parser.add_argument('--salt', type=str, required=True, help='Name of the salt material to use')
args = parser.parse_args()
salt_name = args.salt

# Define ranges for simulation
breeder_thicknesses = np.linspace(0, 50.0, 10)
enrichments = np.linspace(0, 100, 10)

# Initialize results array
results = []

##################################################################################################
# Main Simulation
##################################################################################################
# Loops through all the enrichments for a material for every thickness of the material
for width in breeder_thicknesses:
    for enrich in enrichments:

        ##################################################################################################
        # Materials
        ##################################################################################################

        # Defines the natural berillium reflector
        Breeder = openmc.Material(name="Breeding Layer")
        Breeder.set_density("g/cm3", 1.845)
        Breeder.add_element("Be", 100)
        materials = openmc.Materials([Breeder])

        # Creates though all the salts you want to define below at a constant density of 1.7 g/cm^3

        if salt_name == 'FLiBe':
            salt = openmc.Material(name='FLiBe')
            salt.set_density('g/cm3', 1.7)
            salt.add_element('F', 57.3)
            salt.add_element('Li', 28.2, enrichment=enrich, enrichment_target='Li6')
            salt.add_element('Be', 14.5)

        elif salt_name == 'LiBaBi':
            salt = openmc.Material(name="LiBaBi")
            salt.set_density("g/cm3", 1.7)
            salt.add_element('Li', 20, enrichment=enrich, enrichment_target='Li6')
            salt.add_element('Ba', 10)
            salt.add_element('Bi', 70)

        elif salt_name == 'LiPbBa':
            salt = openmc.Material(name="LiPbBa")
            salt.set_density("g/cm3", 1.7)
            salt.add_element('Li', 25, enrichment=enrich, enrichment_target='Li6')
            salt.add_element('Pb', 60)
            salt.add_element('Ba', 15)

        elif salt_name == 'LiSnZn':
            salt = openmc.Material(name="LiSnZn")
            salt.set_density("g/cm3", 1.7)
            salt.add_element('Li', 65, enrichment=enrich, enrichment_target='Li6')
            salt.add_element('Sn', 25)
            salt.add_element('Zn', 10)

        elif salt_name == 'LiCuPb':
            salt = openmc.Material(name="LiCuPb")
            salt.set_density("g/cm3", 1.7)
            salt.add_element('Li', 40, enrichment=enrich, enrichment_target='Li6')
            salt.add_element('Cu', 20)
            salt.add_element('Pb', 40)

        elif salt_name == 'LiGaPb':
            salt = openmc.Material(name="LiGaPb")
            salt.set_density("g/cm3", 1.7)
            salt.add_element('Li', 35, enrichment=enrich, enrichment_target='Li6')
            salt.add_element('Ga', 10)
            salt.add_element('Pb', 55)

        elif salt_name == 'LiSrPb':
            salt = openmc.Material(name="LiSrPb")
            salt.set_density("g/cm3", 1.7)
            salt.add_element('Li', 30, enrichment=enrich, enrichment_target='Li6')
            salt.add_element('Sr', 50)
            salt.add_element('Pb', 20)

        elif salt_name == 'LiPbZn':
            salt = openmc.Material(name="LiPbZn")
            salt.set_density("g/cm3", 1.7)
            salt.add_element('Li', 30, enrichment=enrich, enrichment_target='Li6')
            salt.add_element('Pb', 60)
            salt.add_element('Zn', 10)

        elif salt_name == 'LiNaSn':
            salt = openmc.Material(name="LiNaSn")
            salt.set_density("g/cm3", 1.7)
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
        settings.batches = 10
        settings.inactive = 0
        settings.particles = 100
        settings.run_mode = 'fixed source'
        settings.source = source
        settings.export_to_xml()
            
        ##################################################################################################
        # Tallies
        ##################################################################################################

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
        with openmc.StatePoint("statepoint.10.h5") as statepoint:
            TBR_mean = statepoint.get_tally(name="TBR").mean.item()
            TBR_std_dev = statepoint.get_tally(name="TBR").std_dev.item()

            # Save results to list
            results.append([width, enrich, TBR_mean, TBR_std_dev])

##################################################################################################
# Postprocessing
##################################################################################################

# Convert results into a NumPy array for structured processing and plotting
results = np.array(results)

# Extract width, enrichment, and TBR values which are needed to generate plots
widths = results[:, 0]
enrichments = results[:, 1]
TBR_values = results[:, 2]

# Reshape arrays for 3D plotting
X, Y = np.meshgrid(np.unique(enrichments), np.unique(widths))
Z = results[:, 2].reshape(len(np.unique(widths)), len(np.unique(enrichments)))

# Identify maximum TBR and the parameter values at this maxima
max_index = np.argmax(results[:, 2])
optimal_width = results[max_index, 0]
optimal_enrichment = results[max_index, 1]
max_tbr = results[max_index, 2]

##################################################################################################
# Plots
##################################################################################################

# Heatmap Plot
fig1, ax1 = plt.subplots(figsize=(10, 7))
heatmap = ax1.contourf(X, Y, Z, levels=20, cmap="viridis")
cbar = plt.colorbar(heatmap)
ax1.set_xlabel("Lithium Enrichment (%)")
ax1.set_ylabel("Breeder Blanket Width (cm)")
ax1.set_title("TBR vs. Width vs. Enrichment (Heatmap)")
cbar.set_label("TBR Value")
plt.savefig(f"tbr_heatmap_{salt_name}.png", dpi=300)
plt.show()

# 3D Plot
fig2 = plt.figure(figsize=(10, 7))
ax2 = fig2.add_subplot(111, projection='3d')
ax2.plot_surface(X, Y, Z, cmap="viridis")
ax2.set_xlabel("Lithium Enrichment (%)")
ax2.set_ylabel("Breeder Blanket Width (cm)")
ax2.set_zlabel("TBR")
ax2.set_title("TBR vs. Width vs. Enrichment (3D Surface Plot)")
plt.savefig(f"tbr_3d_plot{salt_name}.png", dpi=300)
plt.show()

# Plot TBR vs Breeder Thickness at Optimal Enrichment
mask_enrich = results[:, 1] == optimal_enrichment
widths_at_opt_enrich = results[mask_enrich, 0]
tbr_at_opt_enrich = results[mask_enrich, 2]

# Plot TBR vs Breeder Thickness
fig3, ax3 = plt.subplots(figsize=(10, 6))
ax3.plot(widths_at_opt_enrich, tbr_at_opt_enrich, linestyle='-')
ax3.set_xlabel("Breeder Blanket Width (cm)")
ax3.set_ylabel("TBR")
ax3.set_title(f"TBR vs Width at {optimal_enrichment:.2f}% Enrichment for {salt_name}")
plt.grid(True)
plt.savefig(f"tbr_vs_width_at_enrich_{salt_name}.png", dpi=300)
plt.close(fig3)

# Plot TBR vs Enrichment at Optimal Width
mask_width = results[:, 0] == optimal_width
enrichments_at_opt_width = results[mask_width, 1]
tbr_at_opt_width = results[mask_width, 2]

# Plot TBR vs Enrichment
fig4, ax4 = plt.subplots(figsize=(10, 6))
ax4.plot(enrichments_at_opt_width, tbr_at_opt_width, linestyle='--')
ax4.set_xlabel("Lithium Enrichment (%)")
ax4.set_ylabel("TBR")
ax4.set_title(f"TBR vs Enrichment at {optimal_width:.2f} cm Width for {salt_name}")
plt.grid(True)
plt.savefig(f"tbr_vs_enrichment_at_width_{salt_name}.png", dpi=300)
plt.close(fig4)

##################################################################################################
# Save Data
##################################################################################################

# Print to console
print(f"Maximum TBR for {salt_name}: {max_tbr:.5f}")
print(f"Optimal breeder width: {optimal_width:.2f} cm")
print(f"Optimal lithium enrichment: {optimal_enrichment:.2f} %")

# Remove all previosly named text files to prevent rewriting in previous files
for file in glob.glob("raw_data_*.txt"):
    os.remove(file)

# Print results in a readable format
with open(f'raw_data_{salt_name}.txt', 'a') as file:
    file.write("Width (cm) | Enrichment (%) | TBR Mean | TBR Std Dev\n")
    for width, enrichment, tbr_mean, tbr_std in results:
        file.write(f"{width:10.2f} | {enrichment:14.2f} | {tbr_mean:8.5f} | {tbr_std:10.5f}\n")
    file.write(f"\n\n The maximum enrichment for {salt_name} is {np.max(TBR_values)}")
    file.write(f"\n\n The mean enrichment for {salt_name} is {np.mean(TBR_values)}")
    file.write(f"\n\n The median enrichment for {salt_name} is {np.median(TBR_values)}")
    file.write("\n\nMaximum TBR Results:\n")
    file.write(f"Maximum TBR for {salt_name}: {max_tbr:.5f}\n")
    file.write(f"Optimal breeder width: {optimal_width:.2f} m\n")
    file.write(f"Optimal lithium enrichment: {optimal_enrichment:.2f} %\n")

# Remove unneeded files to reduce clutter
os.remove('geometry.xml')
os.remove('materials.xml')
os.remove('settings.xml')
os.remove(f'statepoint.{settings.batches}.h5')
os.remove('summary.h5')
os.remove('tallies.xml')
