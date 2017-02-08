#!/usr/bin/python
import os

### Step 1 creating the required filenames
### NB calling script on Apocrita via ssh requires full pathnames
### NB all path names containing bt**** MUST be changed to relect your individual username

database = "/data/scratch/bt16003/group_project/H_sapiens.db.txt"  # location of the blast database
results = "/data/scratch/bt16003/group_project/results/b*"        # location of the blast outputs
gene_list = "/data/scratch/bt16003/group_project/results/1_gene_list.txt" # name/location of blast gene list outut

uploaded_files = os.listdir('/data/scratch/bt16003/group_project/raw') # names of files unploaded by user
raw_files = []
blast_files = []
blast_results = []


for uploaded_file in uploaded_files:
    # generating full pathnames for raw data
    raw_file = ('/data/scratch/bt16003/group_project/raw/', uploaded_file)
    raw_file = ''.join(raw_file)
    raw_files.append(raw_file)
    # generating full pathnames for primary blast tabular output
    blast_file = ('/data/scratch/bt16003/group_project/temp/blast_', uploaded_file)
    blast_file = ''.join(blast_file)
    blast_files.append(blast_file)
    # generating full pathnames for reformated blast tabular output
    blast_result = ('/data/scratch/bt16003/group_project/results/blast_', uploaded_file)
    blast_result = ''.join(blast_result)
    blast_results.append(blast_result)


### Step 2 running the local blast
for (raw_file, blast_file) in zip(raw_files, blast_files):
    cmd="blastn -db ", database, " -query ", raw_file, \
    \
    " -outfmt \"6 qseqid sseqid stitle\" -max_target_seqs 1 -evalue 0.001 -num_threads 4 -out ", blast_file

    cmd=''.join(cmd)
    print cmd   
    os.system(cmd)
    print (blast_file, "complete")

### Step 3 reformating og blast tabular output
#i) seperate out FPKM value from qseqid
#ii) join stitle components to feature as one column
#iii) Creating an awk array that sums FPKM for redundant hits to the same reference ID (ie isoforms)
for (blast_file, blast_result) in zip(blast_files, blast_results):
    cmd = "cat ", blast_file , " | awk -F ';' '{print $2; next}{print}' | ", \
    \
    "awk '{out= $2 \"\\t\"; for(i = 6; i <= NF; i++){out=out $i\"\\_\"};", \
    \
    "{arr[out] += $1}} END {for (a in arr) {print  a \"\\t\" arr[a]}}' > ", blast_result
    
    cmd= ''.join(cmd)
    os.system(cmd)
    print (blast_result, "complete")


### Step 4 generating a gene reference document for later use in R
#i) concatenate all blast outputs
#ii) sort table such that each reference ID is only included once
#iii) extract the last element between brackets -> this is the gene symbol ie KRAS

cmd= "cat ", results, " | sort -u -k1,1 | awk '{out= $1 \"\\t\" $2; print out}' |", \
\
 " awk -F '[()]' '{out= $1$2$3$4$5 \"\\t\" $2 \"\\t\" $4 ; print out}' | ", \
\
 " awk '{print $1 \"\\t\" $2 \"\\t\" $NF}' > ", gene_list
cmd= ''.join(cmd)
os.system(cmd)


