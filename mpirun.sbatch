#!/bin/bash
#SBATCH --job-name=flamingo_{{RUN_NAME}}
#SBATCH --partition=normal
#SBATCH --output=output-%j.txt
#SBATCH --error=error-%j.txt
#SBATCH --time=48:00:00
#SBATCH -N 1 
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --exclude=grammar082
#SBATCH --exclusive
#SBATCH --mem=0
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
# Your commands go here
module purge
module load gnu14/14.2.0
module load openmpi5/5.0.7
module load hwloc/2.12.0
module load libfabric/1.18.0
module load cmake/3.31.10 
module load gsl/2.8
module load hdf5/1.14.6
module unload ucx
module load ucx/1.18.0-mt
module list

export LIBRARY_PATH=$LIBRARY_PATH:/home/dbi224/local/omp5/fftw3/lib
export LIBRARY_PATH=$LIBRARY_PATH:/home/dbi224/local/omp5/GKlib/lib
export LIBRARY_PATH=$LIBRARY_PATH:/home/dbi224/local/omp5/METIS/lib
export LIBRARY_PATH=$LIBRARY_PATH:/home/dbi224/local/omp5/ParMETIS/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/dbi224/local/omp5/fftw3/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/dbi224/local/omp5/GKlib/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/dbi224/local/omp5/METIS/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/dbi224/local/omp5/ParMETIS/lib

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
export OMP_STACKSIZE=2G
export OMPI_MCA_pml=ucx
export OMPI_MCA_osc=ucx
export UCX_TLS=rc,sm,self
export UCX_RNDV_THRESH=inf
#export OMPI_MCA_pml=cm
#export OMPI_MCA_mtl=ofi
#export FI_PROVIDER='verbs;ofi_rxm'
#export UCX_NET_DEVICES=mlx5_0:1
#export UCX_RNDV_THRESH=inf
#export OMPI_MCA_pml=ob1
#export OMPI_MCA_btl=self,tcp

date


# core dump
ulimit -c 0
ulimit -s unlimited
pwd
mpirun -np 1 --map-by ppr:1:node:PE=32 \
    -x OMP_NUM_THREADS -x OMP_PROC_BIND=false \
    -x OMP_STACKSIZE=2G -x LD_LIBRARY_PATH \
    -x OMPI_MCA_pml -x OMPI_MCA_osc \
    -x UCX_TLS \
    ./swift_mpi  --threads=32 --cosmology --flamingo  --verbose=1   flamingo.yml

