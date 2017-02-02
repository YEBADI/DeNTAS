#!/usr/bin/python

#### NB this script is hosted and executed on apocrita
#### it has been included here as a reference only

import os

### Step 1 creating file names for intermediates

raw_files = os.listdir('/data/scratch/bt16003/group_project/raw')

# creating the file names for intermediates
blast_files =[]
ints_1 = []
ints_2 = []
blast_results = []

for raw_file in raw_files:
	blast_file = ('blast_', raw_file)
	blast_file = ''.join(blast_file)
	blast_files.append(blast_file)

	int_1= ('int_1_blast_', raw_file)
        int_1= ''.join(int_1)
        ints_1.append(int_1)

        int_2= ('int_2_blast_', raw_file)
        int_2= ''.join(int_2)
        ints_2.append(int_2) 

        blast_result = ('blast_', raw_file)
	blast_result = ''.join(blast_result)
	blast_results.append(blast_result)
     

### Step 1 Running local blast
print(" now running blast")

for (raw_file, blast_file) in zip(raw_files, blast_files):
    cmd="blastn -db /data/scratch/bt16003/group_project/reference.fasta.txt -query /data/scratch/bt16003/group_project/raw/", raw_file, \
     " -outfmt \"6 qseqid sseqid evalue stitle\" -max_target_seqs 1 -out /data/scratch/bt16003/group_project/temp/", blast_file
    cmd=''.join(cmd)
    os.system(cmd)
    print (blast_file, "complete")

### this option for Blast (-outfmt "6 qseqid sseqid evalue stitle") will output a table with the following columns
# 1) Query ID, 2) Match ID, 3) E-value, 4) Match description



### Step 2 reformating the blast output table 
## this next piece of code use awk to extract the FPKM value and rearranges the columns into a more user firednly format
# 1) match ID, 2) E-value ID, 3) FPKM, 4) Match description

print("fomatting data output")
  

## this code extracts the FPKM value
for (blast_file , int_1)  in zip(blast_files , ints_1):
    cmd = "awk -F ';' '/^asmbl/{print $2; next}{print}' < /data/scratch/bt16003/group_project/temp/", blast_file, \
     " > /data/scratch/bt16003/group_project/temp/", int_1
    cmd = ''.join(cmd)
    #print cmd
    os.system(cmd)
    #print("done")

## this code reformats the tables
for (int_1 , int_2)  in zip(ints_1, ints_2):
    cmd = "awk \'{out= $2 \"\\t\" $3 \"\\t\" $1 \"\\t\"; for(i=7;i<=NF;i++){out=out" "$i}; print out}' ", \
    " < /data/scratch/bt16003/group_project/temp/" , int_1, " > /data/scratch/bt16003/group_project/temp/", int_2
    cmd = ''.join(cmd)
    os.system(cmd)

print "done"  

#####

print("removing duplicates")

for (int_2 , blast_result) in zip(ints_2 , blast_results):
    cmd= "sort -u -k1,1 /data/scratch/bt16003/group_project/temp/", int_2, " > /data/scratch/bt16003/group_project/results/", blast_result
    cmd= ''.join(cmd)
    os.system(cmd)
    print(blast_result, "done")


##### this code creates the gene list

cmd= "cat /data/scratch/bt16003/group_project/results/b* | sort -u -k1,1 | awk '{out= $1 \"\\t\" $4; print out}' |", \
 " awk -F '[()]' '{out= $1$2$3$4$5 \"\\t\" $2 \"\\t\" $4 ; print out}' | ", \
 " awk '{print $1 \"\\t\" $2 \"\\t\" $NF}' > /data/scratch/bt16003/group_project/results/1_gene_list.txt"
cmd= ''.join(cmd)
os.system(cmd)

### this clears the temp directory of all file intermediates
cmd= 'rm /data/scratch/bt16003/group_project/temp/*'      
os.system(cmd)
