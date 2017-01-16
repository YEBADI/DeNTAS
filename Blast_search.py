#!/usr/bin/python

### Biopython Blast search - progress so far
### - faster/easier to do a local BLAST search
### Step 1 - install blast
### Step 2 - Search NCBI ref_seq RNA by organism for P_Alecto
### Step 3 - Select option to download search results as a fasta file
### Step 4 - make a local blast database

import os

cmd="makeblastdb -in P_Alecto.fasta.txt -parse_seqids -dbtype nucl"
os.system(cmd)

### Step 5 - Data preperations
### unfortunately our sampel data lack the '>' that should be at the start of Fasta ID's therfore not recognized by blast
### Also the files are very large to be used for prototyping

files = ('10s0r1.fasta.txt','12s0r3.fasta.txt', '14s8r2.fasta.txt', \
	'16s24r1.fasta.txt', '18s24r3.fasta.txt', '11s0r2.fasta.txt',\
	'13s8r1.fasta.txt','15s8r3.fasta.txt','17s24r2.fasta.txt')               # these are the unmodified files
results = ('10s0r1.short1.txt','12s0r3.short1.txt', '14s8r2.short1.txt', \
	'16s24r1.short1.txt', '18s24r3.short1.txt', '11s0r2.short1.txt', \
	'13s8r1.short1.txt','15s8r3.short1.txt','17s24r2.short1.txt')            # intermediate files
results_short = ('10s0r1.short.txt','12s0r3.short.txt', '14s8r2.short.txt', \
	'16s24r1.short.txt', '18s24r3.short.txt', '11s0r2.short.txt', \
	'13s8r1.short.txt','15s8r3.short.txt','17s24r2.short.txt')               # shortened file with '>' added to the FASTA ID's


for (file, result) in zip(files, results):
	cmd = "awk '/^asmbl/{print \">\" $1; next}{print}' <", file, ">", result  # adding the >
	cmd = ''.join(cmd)
	os.system(cmd)
	print (file, "done")


for (result, result_short) in zip(results, results_short):	
	cmd = "sed -n '1,1000 p'   ", result, ">", result_short                  # shortening the files
	cmd = ''.join(cmd)
	os.system(cmd)
	print (result, "done")

cmd = "rm 1*short1*"   
os.system(cmd)

### This script should change the Fasta ID to assembled_transcript_1/2/3 etc
### it will also shorten the files to the first 1000 lines (we can adjust this as we like)


### Step 6 running the local Blast alignment


files = ('10s0r1.short.txt','12s0r3.short.txt', '14s8r2.short.txt', \
	'16s24r1.short.txt', '18s24r3.short.txt', '11s0r2.short.txt', \
	'13s8r1.short.txt','15s8r3.short.txt','17s24r2.short.txt')    

results = ('10s0r1.blast.txt','12s0r3.blast.txt', '14s8r2.blast.txt', \
	'16s24r1.blast.txt', '18s24r3.blast.txt', '11s0r2.blast.txt', \
	'13s8r1.blast.txt','15s8r3.blast.txt','17s24r2.blast.txt')    


for (file, result) in zip(files, results):
	cmd="blastn -db reference.fasta.txt -query ", file, " -outfmt \"6 qseqid sseqid evalue stitle\" -max_target_seqs 1 -out results/", result
	print(cmd)
	cmd=''.join(cmd)
	print (result, "complete")

### this option for Blast will output a table with the following columns
# 1) Query ID, 2) Match ID, 3) E-value, 4) Match description

### this next pice of code extracts the FPKM value and rearranges the columns into a more user firednly format
# 1) match ID, 2) E-value ID, 3) FPKM, 4) brief description

files = ('10s0r1.blast.txt' ,'12s0r3.blast.txt' , '14s8r2.blast.txt', \
	'16s24r1.blast.txt', '18s24r3.blast.txt', '11s0r2.blast.txt', \
	'13s8r1.blast.txt','15s8r3.blast.txt','17s24r2.blast.txt')  
intermediates = ('10s0r1.int.txt','12s0r3.int.txt', '14s8r2.int.txt', \
	'16s24r1.int.txt', '18s24r3.int.txt', '11s0r2.int.txt', \
	'13s8r1.int.txt','15s8r3.int.txt','17s24r2.int.txt')  
results = ('10s0r1.txt', '12s0r3.txt','12s0r3.txt', '14s8r2.txt', \
	'16s24r1.txt', '18s24r3.txt', '11s0r2.txt', \
	'13s8r1.txt','15s8r3.txt','17s24r2.txt')   


for (file,intermediate)  in zip(files,intermediates):
	cmd = "awk -F ';' '/^asmbl/{print $2; next}{print}' < ", file, " > ", intermediate
	cmd = ''.join(cmd)
	#print cmd
	os.system(cmd)
    #print("done")

for (intermediate,result)  in zip(intermediates,results):
    cmd = "awk \'{out= $2 \"\\t\" $3 \"\\t\" $1 \"\\t\"; for(i=7;i<=NF;i++){out=out" "$i}; print out}' ", " < " , intermediate, " > ", result
    cmd = ''.join(cmd)
    os.system(cmd)

print "done"  



### the following issues remain with this approach:
### 2) not choosing our own thresholding values
### 3) this approach returns multiple hits for some of the queries (ie transcript variants) 
### 4) I'm having trouble running blast on my computer so have been preparing the data files locally 
###    and then running the blast on apocrita






