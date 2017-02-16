#!/usr/bin/python
import subprocess
import os

### This function can be used to add a ">" to the fasta titles in a mulitfasta if incorrectly formatted
def data_prep(raw_files):
    print("preparing data")
    for raw_file in raw_files:
        cmd = "awk '/^asmbl/{print \">\" $1; next}{print}' < uploads/", raw_file, "> temp/", raw_file  # adding the >
        cmd = ''.join(cmd)
        os.system(cmd)
        print (raw_file, "done")

### The following functions enable the use of Apocrita to run Blast on the users uploaded files remotely

### Copy the contents from DeNTAS' uploads folder to Apocrita's raw folder
def apocrita_upload():
    cmd = "scp ./uploads/* bt16003@login.hpc.qmul.ac.uk:/data/scratch/bt16003/group_project/raw"
    print cmd 
    os.system(cmd)

### The user-selcted organism determins which set of blast scripts are called over ssh
### <organism>_call_blast.sh runs the qsub script <organism>_run _blast.sh which then 
### submits the <organism>_blast.py job to the Apocrita quing system
def apocrita_blast(organism):
    if organism == 'P. alecto':
        cmd = "ssh bt16003@login.hpc.qmul.ac.uk /data/home/bt16003/P_alecto_call_blast.sh"
        print cmd 
        os.system(cmd)
    elif organism == 'H. sapiens':
        cmd = "ssh bt16003@login.hpc.qmul.ac.uk /data/home/bt16003/H_sapiens_call_blast.sh"
        print cmd 
        os.system(cmd)
    elif organism == 'M. musculus':
        cmd = "ssh bt16003@login.hpc.qmul.ac.uk /data/home/bt16003/M_musculus_call_blast.sh"
        print cmd 
        os.system(cmd)

### Copy the contents from Apocrita's results folder to DeNTAS' blast folder
def apocrita_download():
    cmd = "scp bt16003@login.hpc.qmul.ac.uk:/data/scratch/bt16003/group_project/results/* ./blast/"
    os.system(cmd)
    print "blast results downloaded"

### Use ghostscript (a 3rd part suite of pdf manipulation software) to join all
### of the data anlaysis pdf output
def pdf_combine():
    cmd= "gs -q -sPAPERSIZE=a4 -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=static/results/graphs/all_graphs.pdf static/results/graphs/*.pdf"
    os.system(cmd)

def R_analysis(groups):
    x = os.listdir('./blast')[1:len(os.listdir('./blast'))]    # subset to ignore the hidden .DStore file
    files=[]
    for file in x:              # complete path of filenames read into R
        file= ('blast/',file)
        file= ''.join(file)
        files.append(file)
    print x                     # check correct filenames in terminal 
    command = 'Rscript'
    path2script = 'scripts/analysis.R'
    args = files + groups
    cmd = [command, path2script] + args           # constructing the subprocess command 
    subprocess.call(cmd, universal_newlines=True) # begining analysis and parsing in required Args  
    return "analysis complete"

### emptpy the apps local uploads & blast folders and Apocrita's remote raw & results folders
def clean_up():
    cmd = 'rm blast/*' 
    os.system(cmd)
    cmd = 'rm uploads/*' 
    os.system(cmd)
    cmd = "ssh bt16003@login.hpc.qmul.ac.uk 'rm /data/scratch/bt16003/group_project/raw/*' "
    os.system(cmd)
    cmd = "ssh bt16003@login.hpc.qmul.ac.uk 'rm /data/scratch/bt16003/group_project/temp/*' "
    os.system(cmd)
    cmd = "ssh bt16003@login.hpc.qmul.ac.uk 'rm /data/scratch/bt16003/group_project/results/*' "
    os.system(cmd)

### emptpy the local results folders filled by previous runs
def clear_previous_results():
    cmd = 'rm static/results/graphs/*' 
    os.system(cmd)
    cmd = 'rm static/results/tables/*'
    os.system(cmd)
    