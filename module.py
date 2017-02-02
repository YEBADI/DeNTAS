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

def apocrita_upload():
    cmd = "scp ./uploads/*.txt bt16003@login.hpc.qmul.ac.uk:/data/scratch/bt16003/group_project/raw"
    print cmd 
    os.system(cmd)

def apocrita_blast():
    cmd = "ssh bt16003@login.hpc.qmul.ac.uk /data/home/bt16003/call_blast.sh"
    print cmd 
    os.system(cmd)

def apocrita_download():
    cmd = "scp bt16003@login.hpc.qmul.ac.uk:/data/scratch/bt16003/group_project/results/* ./blast/"
    os.system(cmd)
    print "blast results downloaded"

def pdf_combine():
    cmd= "gs -q -sPAPERSIZE=a4 -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=static/results/graphs/all_graphs.pdf static/results/graphs/*.pdf"
    os.system(cmd)

def R_analysis(groups):
    path='./blast'
    x = os.listdir(path)[1:len(os.listdir(path))]
    files=[]
    for file in x:
        file= ('blast/',file)
        file= ''.join(file)
        files.append(file)
    print x
    command = 'Rscript'
    path2script = 'scripts/analysis.R'
    args = files + groups
    cmd = [command, path2script] + args
    subprocess.call(cmd, universal_newlines=True)
    return "analysis complete"

def clean_up():
    cmd = 'rm blast/*' 
    os.system(cmd)
    cmd = 'rm uploads/*' 
    os.system(cmd)
    cmd = "ssh bt16003@login.hpc.qmul.ac.uk 'rm /data/scratch/bt16003/group_project/raw/*' "
    os.system(cmd)
    cmd = "ssh bt16003@login.hpc.qmul.ac.uk 'rm /data/scratch/bt16003/group_project/results/*' "
    os.system(cmd)

