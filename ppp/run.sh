#!/bin/bash
echo "Bash: Running Garfield++"

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <process_index> <dv>"
    exit 1
fi

# Get the process index from the argument
process_index=$1
dv=$2

source /cvmfs/sft.cern.ch/lcg/views/LCG_104c/x86_64-el9-gcc13-opt/setup.sh
source /afs/cern.ch/work/t/tszwarce/garfieldpp/install/share/Garfield/setupGarfield.sh

cd /afs/cern.ch/work/t/tszwarce/ppp/build

./pp $process_index $dv
