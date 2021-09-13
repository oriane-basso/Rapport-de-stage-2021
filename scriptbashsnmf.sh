#!/bin/bash
###Configuration du slurm####
#SBATCH --job-name=Rscript1708smnfpackageLEA
#SBATCH -c 12
#SBATCH --mail-user=oriane.basso@hotmail.com
#SBATCH --mail-type=ALL
#SBATCH -p short
#SBATCH --nodelist=node20
###Par rapport a R
#SBATCH -J lauchRscript
#SBATCH -o output.out
########

#creer un rep perso dans /scratch
cd /scratch
mkdir basso1708-$SLURM_JOB_ID
cd basso1708-$SLURM_JOB_ID

#copie les donnees depuis serveur NAS
scp -r nas:/home/basso/1708falciCHR1 /scratch/basso1708-$SLURM_JOB_ID



#Purge all previously loaded modules
module purge

#load the R module
module load bioinfo/R/4.0.2

#Lancer analyse en se placant au bon endroit
cd 1708falciCHR1

Rscript CHR01scriptR1708.R

#Récupérer les résultats
cd ..
scp -r /scratch/basso1708-$SLURM_JOB_ID/ nas:/home/basso/ResultsmnfLEA1708falciparum

#supprimer mon repertoire perso
cd /scratch
rm -rf basso1708-$SLURM_JOB_ID