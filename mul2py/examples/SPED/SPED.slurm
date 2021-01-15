#!/bin/bash

#SBATCH --partition=GPUQ
#SBATCH --time=00-01:00:00
#SBATCH --job-name="SPED"
#SBATCH --output=SPED-%A.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --mem=64000
#SBATCH --account=share-nv-fys-tem

echo "we are running from this directory: $SLURM_SUBMIT_DIR"
echo "The name of the job is: $SLURM_JOB_NAME"

echo "The job ID is $SLURM_JOB_ID"
echo "The job was run on these nodes: $SLURM_JOB_NODELIST"
echo "Number of nodes: $SLURM_JOB_NUM_NODES"
echo "We are using $SLURM_CPUS_ON_NODE cores"

echo "We are using $SLURM_CPUS_ON_NODE cores per node"
echo "Total of $SLURM_NTASKS cores"

module load foss/2016a
module load CUDA/8.0.61
module load MATLAB/2017a

echo "Running SPED simulation"
matlab -nodisplay -nodesktop -nosplash -r "SPED"

echo "Converting results to HyperSpy format using mul2py"
module load GCCcore/.8.2.0 Python/3.7.2
source /lustre1/projects/itea_lille-nv-fys-tem/MULTEM/mul2py-env/bin/activate
python /lustre1/projects/itea_lille-nv-fys-tem/MULTEM/mul2py/mul2py/examples/convert_ecmat.py SPED_results.ecmat

scontrol show job ${SLURM_JOB_ID} -d