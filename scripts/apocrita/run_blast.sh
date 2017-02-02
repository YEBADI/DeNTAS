#!/bin/sh

#### NB this script is hosted and executed on apocrita
#### it has been included here as a reference only


#$ -cwd                 # use current working directory
#$ -V                   # env variables
#$ -j y # and put all stdout and stderr into the same file 
#$ -l h_rt=1:0:0 # Request 1 hour runtime (upto 240 hours)
#$ -l h_vmem=3G         # Request 1GB RAM / core,

module load blast+
/data/scratch/bt16003/group_project/blast.py
echo 'WORK WORK'  # prints to your output file
