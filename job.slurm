#!/bin/bash
#SBATCH --job-name=openmc_individual
#SBATCH --output=individual_output.log
#SBATCH --error=individual_error.log
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# All Salt Names: FLiBe, LiBaBi, LiPbBa, LiSnZn, LiCuPb, LiGaPb, LiSrPb, LiPbZn, LiNaSn

{
    echo "Running the HYLIFE simulation..."

    cd $SCRATCH/Test_Spheres
    module load OpenMC/0.15.0
    module load MOAB/5.5.1
    module load DAGMC/3.2.3

    salts=("FLiBe" "LiBaBi" "LiPbBa" "LiSnZn" "LiCuPb" "LiGaPb" "LiSrPb" "LiPbZn" "LiNaSn")

    for salt in "${salts[@]}"; do
        echo -e "$salt" | python Infinite_Sphere.py --salt $salt
    done

    echo "Job completed successfully!"
} &>> openmc_output_${SLURM_JOB_ID}.log