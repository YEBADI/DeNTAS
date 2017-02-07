#!/bin/sh
#$ -cwd                 # use current working directory
#$ -V                   # env variables
#$ -j y # and put all stdout and stderr into the same file 
#$ -l h_rt=1:0:0 # Request 1 hour runtime (upto 240 hours)
#$ -l h_vmem=3G         # Request 1GB RAM / core,

# <username> must be changed to reflect your apocrita username eg bt16003
# this script will only work if you have set up a secure key access to apocrita

module load blast+
/data/scratch/<username>/group_project/P_alecto_blast.py
echo 'WORK WORK'  # prints to your output file
