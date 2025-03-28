import openmc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import argparse
import os

# Parse arguments
parser = argparse.ArgumentParser(description='Run OpenMC simulation with specified salt.')
parser.add_argument('--salt', type=str, required=True, help='Name of the salt material to use')
args = parser.parse_args()
salt_name = args.salt

# Define parameter ranges
breeder_thicknesses = np.linspace(0, 50.0, 10)  # Widths from 0m to 10m
enrichments = np.linspace(0, 100, 10)  # Enrichment from 0% to 100% in 10% steps

# Storage for results
results = []

# Loop through Breeder Thicknesses and Enrichments
for width in breeder_thicknesses:
    for enrich in enrichments:

        ##################################################################################################
        # Materials
        ##################################################################################################

        ##################################################################################################
        # Salts
        ##################################################################################################
        
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
            
        # Breeder Material
        Breeder = openmc.Material(name="Breeding Layer")
        Breeder.set_density("g/cm3", 1.845)
        Breeder.add_element("Be", 100)
        materials = openmc.Materials([Breeder])

        materials.append(salt)
        materials.export_to_xml()

        ##################################################################################################
        # Geometry
        ##################################################################################################

        # Set Radii
        inner_radius, outer_radius = 6, 2000
        reflector_radius = 6 + width
        boundary = outer_radius + 1

        # Sphere within a Aphere with no Leakage
        sphere_inner = openmc.Sphere(r=inner_radius)
        sphere_reflector = openmc.Sphere(r=reflector_radius)
        sphere_outer = openmc.Sphere(r=outer_radius)
        sphere_leakage = openmc.Sphere(r=boundary, boundary_type='reflective')

        # Define Cells
        center_void = -sphere_inner
        region_reflector = +sphere_inner & -sphere_reflector
        region_salt = +sphere_reflector & -sphere_outer
        outer_void = +sphere_outer & -sphere_leakage

        ##################################################################################################
        # Cells
        ##################################################################################################

        # Define and Fill Cells
        center_cell = openmc.Cell(fill=None, region=center_void)
        cell_reflector = openmc.Cell(fill=Breeder, region=region_reflector)
        cell_salt = openmc.Cell(fill=salt, region=region_salt)
        outer_cell = openmc.Cell(fill=None, region=outer_void)

        # Define Universe
        universe = openmc.Universe(cells=[center_cell, cell_salt, cell_reflector, outer_cell])
        geometry = openmc.Geometry(universe)
        geometry.export_to_xml()

        ##################################################################################################
        # Settings
        ##################################################################################################

        # Source Definition
        source = openmc.IndependentSource()
        source.particle = 'neutron'
        source.energy = openmc.stats.Discrete([14.1e6], [1.0])
        source.angle = openmc.stats.Isotropic()
        source.space = openmc.stats.Point((0.0, 0.0, 0.0))

        # Settings
        settings = openmc.Settings()
        settings.batches = 10
        settings.inactive = 0
        settings.particles = 10
        settings.run_mode = 'fixed source'
        settings.source = source
        settings.export_to_xml()
            
        ##################################################################################################
        # Tallies
        ##################################################################################################

        tallies = openmc.Tallies()

        # TBR Tally for the Salt being Run
        TBR_tally = openmc.Tally(name='TBR')
        TBR_tally.filters = [openmc.MaterialFilter([salt])]
        TBR_tally.scores = ['(n,Xt)']
        tallies.append(TBR_tally)

        tallies.export_to_xml()

        ##################################################################################################
        # Run OpenMC
        ##################################################################################################
        
        openmc.run()

        ##################################################################################################
        # Extract Results
        ##################################################################################################
        
        # Take the data and append it to the results
        with openmc.StatePoint("statepoint.10.h5") as statepoint:
            TBR_mean = statepoint.get_tally(name="TBR").mean.item()  # Ensure it's a scalar
            TBR_std_dev = statepoint.get_tally(name="TBR").std_dev.item()  # Ensure it's a scalar

            # Save results to list
            results.append([width, enrich, TBR_mean, TBR_std_dev])

# Convert results into a NumPy array
results = np.array(results)

# Extract width, enrichment, and TBR values
widths = results[:, 0]  # Y-axis
enrichments = results[:, 1]  # X-axis
TBR_values = results[:, 2]  # Z-axis (TBR)

# Reshape arrays for 3D plotting
X, Y = np.meshgrid(np.unique(enrichments), np.unique(widths))
Z = results[:, 2].reshape(len(np.unique(widths)), len(np.unique(enrichments)))

# Heatmap Plot
fig1, ax1 = plt.subplots(figsize=(10, 7))
heatmap = ax1.contourf(X, Y, Z, levels=20, cmap="coolwarm")
cbar = plt.colorbar(heatmap)
ax1.set_xlabel("Lithium Enrichment (%)")
ax1.set_ylabel("Breeder Blanket Width (m)")
ax1.set_title("TBR vs. Width vs. Enrichment (Heatmap)")
cbar.set_label("TBR Value")
plt.savefig(f"tbr_heatmap_{salt_name}.png", dpi=300)
plt.show()

# 3D Plot
fig2 = plt.figure(figsize=(10, 7))
ax2 = fig2.add_subplot(111, projection='3d')
ax2.plot_surface(X, Y, Z, cmap="viridis")
ax2.set_xlabel("Lithium Enrichment (%)")
ax2.set_ylabel("Breeder Blanket Width (m)")
ax2.set_zlabel("TBR")
ax2.set_title("TBR vs. Width vs. Enrichment (3D Surface Plot)")
plt.savefig(f"tbr_3d_plot{salt_name}.png", dpi=300)
plt.show()

print("Both heatmap and 3D plot have been saved.")

# Print results in a readable format
with open(f'raw_data_{salt_name}.txt', 'a') as file:
    file.write("Width (m) | Enrichment (%) | TBR Mean | TBR Std Dev\n")
    for width, enrichment, tbr_mean, tbr_std in results:
        file.write(f"{width:8.2f} | {enrichment:12.2f} | {tbr_mean:8.5f} | {tbr_std:10.5f}\n")
    file.write(f"\n\n The maximum enrichment for {salt_name} is {np.max(TBR_values)}")
    file.write(f"\n\n The mean enrichment for {salt_name} is {np.mean(TBR_values)}")
    file.write(f"\n\n The median enrichment for {salt_name} is {np.median(TBR_values)}")

# Remove all unneded files to reduce clutter
os.remove('geometry.xml')
os.remove('materials.xml')
os.remove('settings.xml')
os.remove(f'statepoint.{settings.batches}.h5')
os.remove(f'summary.h5')
os.remove(f'tallies.xml')
