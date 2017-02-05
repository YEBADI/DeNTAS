#### NB this script is hosted and executed on apocrita

#$ -pe smp 4       # Request 4 CPU cores (please add this to the script executed in apocrita so we can get 4 cores)
#set the number GB per core to 1 please

#### it has been included here as a reference only

export SGE_ROOT="/opt/sge"
/opt/sge/bin/lx-amd64/qsub -sync y /data/home/bt16003/run_blast.sh

