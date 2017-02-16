# <username> must be changed to reflect your apocrita username eg bt16003
# this script will only work if you have set up a secure key access to apocrita

export SGE_ROOT="/opt/sge"
/opt/sge/bin/lx-amd64/qsub -sync y /data/home/<username>/M_musculus_run_blast.sh

