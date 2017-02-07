#!/bin/sh

# <username> must be changed to reflect your apocrita username eg bt16003
# this script will only work if you have set up a secure key access to apocrita

### copying over the db FASTAs
scp ./H_sapiens.db.txt.gz <username>@login.hpc.qmul.ac.uk:
scp ./P_alecto.db.txt.gz <username>@login.hpc.qmul.ac.uk:

### copying over the required scripts
scp ./call_blast_setup.sh <username>@login.hpc.qmul.ac.uk:
scp ./blast_setup.sh <username>@login.hpc.qmul.ac.uk:
scp ./H_sapiens_call_blast.sh <username>@login.hpc.qmul.ac.uk:
scp ./H_sapiens_run_blast.sh <username>@login.hpc.qmul.ac.uk:
scp ./P_alecto_call_blast.sh <username>@login.hpc.qmul.ac.uk:
scp ./P_alecto_run_blast.sh <username>@login.hpc.qmul.ac.uk:
scp ./P_alecto_blast.py <username>@login.hpc.qmul.ac.uk:
scp ./H_sapiens_blast.py <username>@login.hpc.qmul.ac.uk:

### running the blast setup
# 1) creates the correct directory structure in your scratch space
# 2) 
ssh <username>@login.hpc.qmul.ac.uk ./call_blast_setup.sh