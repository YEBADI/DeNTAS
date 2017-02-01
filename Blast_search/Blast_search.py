#running blast off of apocrita
#ssh [secret key] @login.hpc.qmul.ac.uk
#scp [upload files]

#!/usr/bin/python

### Biopython Blast search 
### - faster/easier to do a local BLAST search
### Step 1 - install blast
### Step 2 - Search NCBI ref_seq RNA by organism for P_Alecto
### Step 3 - Select option to download search results as a fasta file
### Step 4 - make a local blast database
### NB blast doesn't run on my computer so i have been uploading the raw data to Apocrita....
### ....running this script and then downloading the results 

### to run this script you require
### i) a local blast db in the working directory ii) a subdirectory raw/ containg the unmodified fastas
### ii) a subdirectory results/ where the final, reformated output will appear


import os
### Step 4 making a local blast DB

#cmd="makeblastdb -in reference.fasta.txt -parse_seqids -dbtype nucl"
#os.system(cmd)

### Step 5 - Data preperations
### Our data lacks the '>' that should be at the start of Fasta ID's the following code...
### ...adds it in

files = ('10s0r1.fasta.txt','12s0r3.fasta.txt', '14s8r2.fasta.txt', \
	'16s24r1.fasta.txt', '18s24r3.fasta.txt', '11s0r2.fasta.txt',\
	'13s8r1.fasta.txt','15s8r3.fasta.txt','17s24r2.fasta.txt')               # these are the raw files
results = ('10s0r1.fasta1.txt','12s0r3.fasta1.txt', '14s8r2.fasta1.txt', \
	'16s24r1.fasta1.txt', '18s24r3.fasta1.txt', '11s0r2.fasta1.txt', \
	'13s8r1.fasta1.txt','15s8r3.fasta1.txt','17s24r2.fasta1.txt')            # files with '>' added


for (file, result) in zip(files, results):
	cmd = "awk '/^asmbl/{print \">\" $1; next}{print}' < raw/", file, ">", result  # adding the >
	cmd = ''.join(cmd)
	os.system(cmd)
	print (file, "done")

### Step 6 running the local Blast alignment

print(" now running blast")

files = ('10s0r1.fasta1.txt','12s0r3.fasta1.txt', '14s8r2.fasta1.txt', \
	'16s24r1.fasta1.txt', '18s24r3.fasta1.txt', '11s0r2.fasta1.txt', \
	'13s8r1.fasta1.txt','15s8r3.fasta1.txt','17s24r2.fasta1.txt')    

results = ('10s0r1.blast.txt','12s0r3.blast.txt', '14s8r2.blast.txt', \
	'16s24r1.blast.txt', '18s24r3.blast.txt', '11s0r2.blast.txt', \
	'13s8r1.blast.txt','15s8r3.blast.txt','17s24r2.blast.txt')    


for (file, result) in zip(files, results):
	cmd="blastn -db reference.fasta.txt -query ", file, " -outfmt \"6 qseqid sseqid evalue stitle\" -max_target_seqs 1 -out", result
	cmd=''.join(cmd)
	os.system(cmd)
	print (result, "complete")

# when you call blastn
#blastn -query /path/to/myquery.fasta -db /path/to/db -evalue 0.001 -out /path/to/myoutput.tsv -outfmt "6 qacc sacc sseqid evalue qstart qend sstart send"
	
### this option for Blast (-outfmt "6 qseqid sseqid evalue stitle") will output a table with the following columns
# 1) Query ID, 2) Match ID, 3) E-value, 4) Match description



### Step 6 reformating the blast output table 

### this next piece of code use awk to extract the FPKM value and rearranges the columns into a more user firednly format
# 1) match ID, 2) E-value ID, 3) FPKM, 4) Match description

print("fomatting data output")

files = ('10s0r1.blast.txt' ,'12s0r3.blast.txt' , '14s8r2.blast.txt', \
	'16s24r1.blast.txt', '18s24r3.blast.txt', '11s0r2.blast.txt', \
	'13s8r1.blast.txt','15s8r3.blast.txt','17s24r2.blast.txt')  
intermediates = ('10s0r1.int.txt','12s0r3.int.txt', '14s8r2.int.txt', \
	'16s24r1.int.txt', '18s24r3.int.txt', '11s0r2.int.txt', \
	'13s8r1.int.txt','15s8r3.int.txt','17s24r2.int.txt')  
results = ('10s0r1.full.txt', '12s0r3.full.txt', '14s8r2.full.txt', \
	'16s24r1.full.txt', '18s24r3.full.txt', '11s0r2.full.txt', \
	'13s8r1.full.txt','15s8r3.full.txt','17s24r2.full.txt')   

## this code extracts the FPKM value
for (file,intermediate)  in zip(files,intermediates):
	cmd = "awk -F ';' '/^asmbl/{print $2; next}{print}' < ", file, " > ", intermediate
	cmd = ''.join(cmd)
	#print cmd
	os.system(cmd)
    #print("done")

## this code reformats the tables
for (intermediate,result)  in zip(intermediates,results):
    cmd = "awk \'{out= $2 \"\\t\" $3 \"\\t\" $1 \"\\t\"; for(i=7;i<=NF;i++){out=out" "$i}; print out}' ", " < " , intermediate, " > ", result
    cmd = ''.join(cmd)
    os.system(cmd)

print "done"  

#####

print("removing duplicates")

files = ('10s0r1.full.txt', '12s0r3.full.txt', '14s8r2.full.txt', \
	'16s24r1.full.txt', '18s24r3.full.txt', '11s0r2.full.txt', \
	'13s8r1.full.txt','15s8r3.full.txt','17s24r2.full.txt')   

results = ('10s0r1.txt', '12s0r3.txt', '14s8r2.txt', \
	'16s24r1.txt', '18s24r3.txt', '11s0r2.txt', \
	'13s8r1.txt','15s8r3.txt','17s24r2.txt')


for (file,result) in zip(files,results):
	cmd= "sort -u -k1,1 ", file, " > results/", result
	cmd= ''.join(cmd)
	os.system(cmd)
	print(file, "done")

### this clears the working directory of all file intermediates
cmd= 'rm 1*'      
os.system(cmd)

### Step 8 creating a gene 'look-up' table
### a gene_list identification table must be made to enable looking up of gene by ref_seqID

# the following commands can be cut and past into the linux command line
#making a gene ID - function refrence list
#1) cat 1* | sort -u -k1,1 > functions.txt
#2) awk '{out= $1 "\t" $4; print out}' < functions.txt > functions2.txt
#3) awk -F '[()]' '{out= $1$2$3$4$5 "\t" $2; print out}' < functions2.txt > gene_list.txt
#4) rm  functions*

#or in one line
# cat 1* | sort -u -k1,1 | awk '{out= $1 "\t" $4; print out}' | awk -F '[()]' '{out= $1$2$3$4$5 "\t" $2; print out}' >gene_list.txt




### the following issues remain with this approach:
### 1) not choosing our own thresholding values !!!!
### 2) the blast output table has multiple queries mapping to the same ref_seq ID
### 3) this approach just collapse the ref_seqID so that that column only contains unique values
### ... it would be good to somehow add all of the FPKM values for a single ref_seqID
### 4) this approach does not include any of the transcripts that do not map to a ref_seqID - ...
###... even if they were to be present in all samples - ie maybe the Hendra genes?


