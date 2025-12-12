#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=180G
#SBATCH --mem-per-cpu=1024M
#SBATCH --tmp=24G
#SBATCH --job-name=spectra
#SBATCH --partition=akanane
#SBATCH --time=7-0
#SBATCH --output=job.out
#SBATCH --error=job.err
#SBATCH --mail-user='aarti@udel.edu'
#SBATCH --mail-type=END,FAIL
#SBATCH --export=NONE

cd $SLURM_SUBMIT_DIR

#cp *.inp $SCRATCHDIR
#cp job.out $SCRATCHDIR

#cd $SCRATCHDIR
./job.out

#cp output* $SGE_O_WORKDIR
#cp fort.* $SGE_O_WORKDIR
#cp *.out $SGE_O_WORKDIR
#cd ..
#rm -r $SCRATCHDIR

