#!/bin/bash

# Print usage if the arguments given are not correct.
if [[ $# -eq 0 ]] ; then
    echo ''
    echo 'Usage: RunAllPsiIns [PSI4EXE] [THREADS]'
    echo ''
    echo 'Example: RunAllPsiIns "~/Anaconda3/envs/Psi4/bin/psi4" 2'
    echo ''
    echo 'Uses the Psi4 executable to run all jobs in the current directory.'
    echo ''
    echo 'Mandatory arguments:'
    echo '  [PSI4EXE] Location of the Psi4 executable, e.g. "/usr/bin/psi4". (Default is psi4)'
    echo ''
    echo 'Optional arguments:'
    echo '  [THREADS] Number of execution threads, e.g. 2. (Default is 4)'
    echo ''
    exit 0
fi

# Mandatory argument 1 is the Psi4 executable location
ARG1=${1}

# Optional argument 2 is the number threads to use.
ARG2=${2:-4}

for FILE in *.in
do
    NOEXTNAME=`echo $FILE | cut -d'.' --complement -f2-`
    echo "Executing ${NOEXTNAME}.in with ${ARG1} on ${ARG2} threads."
    #psi4 -i ${NOEXTNAME}.in -o ${NOEXTNAME}.out -n 4
    ${ARG1} -i ${NOEXTNAME}.in -o ${NOEXTNAME}.out -n ${ARG2}
done
