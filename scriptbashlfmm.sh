#!/bin/bash
###Configuration du slurm####
#SBATCH --job-name=Rscript1708smnfpackageLEA
#SBATCH -c 2
#SBATCH --mail-user=oriane.basso@hotmail.com
#SBATCH --mail-type=ALL
#SBATCH -p supermem
###Par rapport a R
#SBATCH -J lauchRscript
#SBATCH -o output.out
########

CHR=1
datedujour=2308
env="latVA"
sp="falci"
datecreationscriptR=2208
Kstart=6
Kend=12
nomvariableenv="latVAB"



#Creer file resultat
cd
mkdir resultsEA${env}${sp}${datedujour}
cd resultsEA${env}${sp}${datedujour}
mkdir CHR${CHR}
cd

#creer un rep perso dans /scratch
cd /scratch
mkdir basso${datedujour}-$SLURM_JOB_ID
cd basso${datedujour}-$SLURM_JOB_ID

#copie les donnees depuis serveur NAS
scp -r nas:/home/basso/falciCHR$CHR /scratch/basso${datedujour}-$SLURM_JOB_ID
cd falciCHR$CHR


#Purge all previously loaded modules
module purge

#load the R module
module load bioinfo/R/4.0.2

#Lancer analyse en se placant au bon endroit
#Boucle pour lancer analyse
for i in {6..12}; do
Rscript SCRIPTCHR${CHR}${env}B${datecreationscriptR}K${i}.R
mkdir K${i}${env}B
mv *lfmm* K${i}${env}B/
mv plotpv.pdf plotpvK${i}.pdf
mv plotpvK${i}.pdf K${i}${env}B/
mv *multiple* K${i}${env}B/
mv *candi* K${i}${env}B/
rm Chr_0${CHR}.removed
mkdir ${datedujour}allK${env}B
mv K${i}${env}B/ ${datedujour}allK${env}B/;done



#Récupérer les résultats
cd ..
scp -r /scratch/basso${datedujour}-$SLURM_JOB_ID/ nas:/home/basso/resultsEA${env}${sp}${datedujour}/CHR${CHR}

#supprimer mon repertoire perso
cd /scratch
rm -rf basso${datedujour}-$SLURM_JOB_ID