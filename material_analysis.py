# NOTES: MATERIALS THAT ARE COMMENTED OUT ARE FOR ACTUAL RUNS. CURRENTLY SILL 
# BEING TESTED SO DO NOT INCREASE BATCHES OR STEPS UNTIL AFTER INITIAL ANALYSIS.
# MATERIALS ALSO NEED TO HAVE THEIR CORRECT DENSITIES SET AND HAVE THE add_s_alpha_beta()
# ADDED FOR RELEVANT FUNCTIONS

import openmc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define parameter ranges
breeder_thicknesses = np.linspace(0, 10.0, 10)  # Widths from 0m to 10m
enrichments = np.linspace(0, 100, 10)  # Enrichment from 0% to 100% in 10% steps

# Storage for results
results = []

for width in breeder_thicknesses:
    for enrichment in enrichments:

        ##################################################################################################
        # Materials
        ##################################################################################################

        Breeder = openmc.Material(name="Breeding Layer")
        Breeder.set_density("g/cm3", 1.845)
        Breeder.add_element("Be", 100)
        materials = openmc.Materials([Breeder])

        FLiBe = openmc.Material(name='FLiBe')
        FLiBe.set_density('g/cm3', 1.937)
        FLiBe.add_element('F', 57.3)
        FLiBe.add_element('Li', 28.2, enrichment=enrichment, enrichment_target='Li6')
        FLiBe.add_element('Be', 14.5)
        materials.append(FLiBe)

        # LiBaBi = openmc.Material(name="LiBaBi")
        # LiBaBi.set_density("g/cm3", 1.7)
        # LiBaBi.add_element('Li', 20, enrichment=enrich, enrichment_target='Li6')
        # LiBaBi.add_element('Ba', 10)
        # LiBaBi.add_element('Bi', 70)
        # testing_materials.append(LiBaBi)

        # LiPbBa = openmc.Material(name="LiPbBa")
        # LiPbBa.set_density("g/cm3", 1.7)
        # LiPbBa.add_element('Li', 25, enrichment=enrich, enrichment_target='Li6')
        # LiPbBa.add_element('Pb', 60)
        # LiPbBa.add_element('Ba', 15)
        # testing_materials.append(LiPbBa)

        # LiSnZn = openmc.Material(name="LiSnZn")
        # LiSnZn.set_density("g/cm3", 1.7)
        # LiSnZn.add_element('Li', 65, enrichment=enrich, enrichment_target='Li6')
        # LiSnZn.add_element('Sn', 25)
        # LiSnZn.add_element('Zn', 10)
        # testing_materials.append(LiSnZn)

        # LiCuPb = openmc.Material(name="LiCuPb")
        # LiCuPb.set_density("g/cm3", 1.7)
        # LiCuPb.add_element('Li', 40, enrichment=enrich, enrichment_target='Li6')
        # LiCuPb.add_element('Cu', 20)
        # LiCuPb.add_element('Pb', 40)
        # testing_materials.append(LiCuPb)

        # LiGaPb = openmc.Material(name="LiGaPb")
        # LiGaPb.set_density("g/cm3", 1.7)
        # LiGaPb.add_element('Li', 35, enrichment=enrich, enrichment_target='Li6')
        # LiGaPb.add_element('Ga', 10)
        # LiGaPb.add_element('Pb', 55)
        # testing_materials.append(LiGaPb)

        # LiSrPb = openmc.Material(name="LiSrPb")
        # LiSrPb.set_density("g/cm3", 1.7)
        # LiSrPb.add_element('Li', 30, enrichment=enrich, enrichment_target='Li6')
        # LiSrPb.add_element('Sr', 50)
        # LiSrPb.add_element('Pb', 20)
        # testing_materials.append(LiSrPb)

        # LiPbZn = openmc.Material(name="LiPbZn")
        # LiPbZn.set_density("g/cm3", 1.7)
        # LiPbZn.add_element('Li', 30, enrichment=enrich, enrichment_target='Li6')
        # LiPbZn.add_element('Pb', 60)
        # LiPbZn.add_element('Zn', 10)
        # testing_materials.append(LiPbZn)

        # LiNaSn = openmc.Material(name="LiNaSn")
        # LiNaSn.set_density("g/cm3", 1.7)
        # LiNaSn.add_element('Li', 55, enrichment=enrich, enrichment_target='Li6')
        # LiNaSn.add_element('Na', 30)
        # LiNaSn.add_element('Sn', 15)
        # testing_materials.append(LiNaSn)

        materials.export_to_xml()

        ##################################################################################################
        # Geometry
        ##################################################################################################

        inner_radius, outer_radius = 6, 2000
        reflector_radius = 6 + width
        boundary = outer_radius + 1

        sphere_inner = openmc.Sphere(r=inner_radius)
        sphere_reflector = openmc.Sphere(r=reflector_radius)
        sphere_outer = openmc.Sphere(r=outer_radius)
        sphere_leakage = openmc.Sphere(r=boundary, boundary_type='reflective')

        center_void = -sphere_inner
        region_reflector = +sphere_inner & -sphere_reflector
        region_salt = +sphere_reflector & -sphere_outer
        outer_void = +sphere_outer & -sphere_leakage

        ##################################################################################################
        # Cells
        ##################################################################################################

        center_cell = openmc.Cell(fill=None, region=center_void)
        cell_reflector = openmc.Cell(fill=Breeder, region=region_reflector)
        cell_salt = openmc.Cell(fill=FLiBe, region=region_salt)
        outer_cell = openmc.Cell(fill=None, region=outer_void)

        universe = openmc.Universe(cells=[center_cell, cell_salt, cell_reflector, outer_cell])
        geometry = openmc.Geometry(universe)

        geometry.export_to_xml()

        ##################################################################################################
        # Settings
        ##################################################################################################

        source = openmc.IndependentSource()
        source.particle = 'neutron'
        source.energy = openmc.stats.Discrete([14.1e6], [1.0])
        source.angle = openmc.stats.Isotropic()
        source.space = openmc.stats.Point((0.0, 0.0, 0.0))

        settings = openmc.Settings()
        settings.batches = 10
        settings.inactive = 10
        settings.particles = 10
        settings.run_mode = 'fixed source'
        settings.source = source
        settings.export_to_xml()
            
        ##################################################################################################
        # Tallies
        ##################################################################################################

        tallies = openmc.Tallies()

        TBR_tally = openmc.Tally(name='TBR')
        TBR_tally.filters = [openmc.MaterialFilter([FLiBe])]
        TBR_tally.scores = ['(n,Xt)']
        tallies.append(TBR_tally)

        tallies.export_to_xml()

        # Run OpenMC simulation
        openmc.run()

        # Extract results
        with openmc.StatePoint("statepoint.10.h5") as statepoint:
            TBR_mean = statepoint.get_tally(name="TBR").mean.item()  # Ensure it's a scalar
            TBR_std_dev = statepoint.get_tally(name="TBR").std_dev.item()  # Ensure it's a scalar

            # Save results to list
            results.append([width, enrichment, TBR_mean, TBR_std_dev])

# Convert results into a NumPy array for structured processing
results = np.array(results)

# Extract width, enrichment, and TBR values
widths = results[:, 0]  # Y-axis
enrichments = results[:, 1]  # X-axis
TBR_values = results[:, 2]  # Z-axis (TBR)

# Reshape arrays for 3D plotting
X, Y = np.meshgrid(np.unique(enrichments), np.unique(widths))
Z = results[:, 2].reshape(len(np.unique(widths)), len(np.unique(enrichments)))

# 3D Surface Plot
fig, ax = plt.subplots(figsize=(10, 7))
heatmap = ax.contourf(X, Y, Z, levels=20, cmap="coolwarm")  # Blue = Low, Red = High
cbar = plt.colorbar(heatmap)

ax.set_xlabel("Lithium Enrichment (%)")
ax.set_ylabel("Breeder Blanket Width (m)")
ax.set_title("TBR vs. Width vs. Enrichment (Heatmap)")
cbar.set_label("TBR Value")
ax.set_title("TBR vs. Width vs. Enrichment")

# Save plot
plt.savefig("tbr_heatmap.png", dpi=300)
plt.show()

print("Plot saved as 'tbr_surface_plot.png'")

# Print results in a readable format
print("Width (m) | Enrichment (%) | TBR Mean | TBR Std Dev")
for width, enrichment, tbr_mean, tbr_std in results:
    print(f"{width:8.2f} | {enrichment:12.2f} | {tbr_mean:8.5f} | {tbr_std:10.5f}")

os.remove('geometry.xml')
os.remove('materials.xml')
os.remove('settings.xml')
os.remove(f'statepoint.{settings.batches}.h5')
os.remove(f'summary.h5')
os.remove(f'tallies.xml')