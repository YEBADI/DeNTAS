#!/usr/bin/python
import subprocess
import os


def data_prep(raw_files):
    print("preparing data")

    for raw_file in raw_files:
        cmd = "awk '/^asmbl/{print \">\" $1; next}{print}' < uploads/", raw_file, "> temp/", raw_file  # adding the >
        cmd = ''.join(cmd)
        os.system(cmd)
        print (raw_file, "done")

def run_blast(raw_files):
    print(" now running blast")
    blast_results=[]
    for file in raw_files:
        result = ('blast_', file)
        result = ''.join(result)
        blast_results.append(result)

    for (file, blast_result) in zip(raw_files, blast_results):
        cmd="blastn -query ", file, " -db reference.fasta.txt -evalue 0.001 -outfmt \"6 qseqid sseqid evalue stitle\" -max_target_seqs 1 -out", result
        cmd=''.join(cmd)
        os.system(cmd)
        print (blast_result, "complete")
    return blast_results

def format_blast(blast_results):
    print("fomatting data output")
    
    ints_1=[]
    for blast_result in blast_results:
        int_1= ('int_1_', blast_result)
        int_1= ''.join(int_1)
        ints_1.append(int_1)  
    
    ints_2=[]
    for blast_result in blast_results:
        int_2= ('int_2_', blast_result)
        int_2= ''.join(int_2)
        ints_2.append(int_2)

    ## this code extracts the FPKM value
    for (blast_result, int_1)  in zip(blast_results, ints_1):
        cmd = "awk -F ';' '/^asmbl/{print $2; next}{print}' < temp/", blast_result, " > temp/", int_1
        cmd = ''.join(cmd)
        os.system(cmd)

    ## this code reformats the tables
    for (int_1, int_2)  in zip( ints_1, ints_2):
        cmd = "awk \'{out= $2 \"\\t\" $3 \"\\t\" $1 \"\\t\"; for(i=7;i<=NF;i++){out=out " "$i \"_\" }; print out}' < temp/" , int_1, " > temp/", int_2
        cmd = ''.join(cmd)
        print (int_2, "done")
        os.system(cmd)

    ## this code removes duplicates & saves output in blast/ for use by R
    print("removing duplicates")

    for (int_2,blast_result) in zip(ints_2, blast_results):
        cmd= "sort -u -k1,1 temp/", int_2, " > blast/", blast_result
        cmd= ''.join(cmd)
        os.system(cmd)
        print(blast_result, "done")

    ### this clears the temp folder of intermediates
    cmd= 'rm temp/int*'      
    os.system(cmd)
    cmd= 'rm temp/1*'      
    os.system(cmd)

def genelist():
    cmd= "cat blast/b* | sort -u -k1,1 | awk '{out= $1 \"\\t\" $4; print out}' | awk -F '[()]' '{out= $1$2$3$4$5 \"\\t\" $2 \"\\t\" $4 ; print out}' | awk '{print $1 \"\\t\" $2 \"\\t\" $NF}' > blast/1_gene_list.txt"
    os.system(cmd)


def R_analysis(groups):
    path='./blast'
    x = os.listdir(path)[1:len(os.listdir(path))]
    files=[]
    for file in x:
        file= ('blast/',file)
        file= ''.join(file)
        files.append(file)

    command = 'Rscript'
    path2script = 'scripts/analysis.R'
    args = files + groups

    cmd = [command, path2script] + args
    subprocess.call(cmd, universal_newlines=True)
    return "analysis complete"


