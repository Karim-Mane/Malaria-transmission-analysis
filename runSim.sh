#!/usr/bin/env bash

set -e;
set -o pipefail;

#echo $1
JID1=`uuidgen`
bsub -q normal -a openmpi -J "$JID1[1-700]%50" -n 2 -M1000 -R"select[mem>1000] rusage[mem=1000]" -oo runSim.out -eo runSim.err -G team273-vwork ./getRelatedness.sh --yaml-file $1 
sleep 2
#JID2=`uuidgen`
#bsub -w "done($JID1)||exit($JID1)" -q normal -J "$JID2" -n 1 -oo attend.out -eo attend.err -G team273-vwork ./Attente.sh 
