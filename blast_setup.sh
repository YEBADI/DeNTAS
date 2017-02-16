#!/bin/sh
#$ -cwd                 # use current working directory
#$ -V                   # env variables
#$ -l h_vmem=1G         # Request 1GB RAM / core,

# <username> must be changed to reflect your apocrita username eg bt16003
# this script will only work if you have set up a secure key access to apocrita

module load blast+

# creatibng the required directories in your scratch space
mkdir /data/scratch/<username>/group_project
mkdir /data/scratch/<username>/group_project/raw
mkdir /data/scratch/<username>/group_project/temp
mkdir /data/scratch/<username>/group_project/results

# Transfering the required scripts and FASTAs to your scratch space
mv P_alecto.db.txt.gz /data/scratch/<username>/group_project
mv H_sapiens.db.txt.gz /data/scratch/<username>/group_project
mv M_musculus.db.txt.gz /data/scratch/<username>/group_project

mv P_alecto_blast.py /data/scratch/<username>/group_project
mv H_sapiens_blast.py /data/scratch/<username>/group_project
mv M_musculus_blast.py /data/scratch/<username>/group_project

# Unziping the compressed FASTAs
gunzip /data/scratch/<username>/group_project/P_alecto.db.txt.gz
gunzip /data/scratch/<username>/group_project/H_sapiens.db.txt.gz
gunzip /data/scratch/<username>/group_project/M_musculus.db.txt.gz

# Creating the local blast databases in your scratch space
makeblastdb -in /data/scratch/<username>/group_project/P_alecto.db.txt -parse_seqids -dbtype nucl
makeblastdb -in /data/scratch/<username>/group_project/H_sapiens.db.txt -parse_seqids -dbtype nucl
makeblastdb -in /data/scratch/<username>/group_project/M_musculus.db.txt -parse_seqids -dbtype nucl

