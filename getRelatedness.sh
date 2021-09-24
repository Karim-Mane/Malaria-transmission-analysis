#!/bin/bash

module load R


function usage 
{
	cat << HEREDOC
	
	 About:	run simulation.
	 
	 Usage:	$progname [--yaml-file YAML_FILE] in interactive mode
	 or bsub -q long -a openmpi -J "runSim[1-700]" -n 20 -M1000 -R"select[mem>1000] rusage[mem=1000]" -oo runSim.out -eo runSim.err -G team273-vwork 
	 ./$progname [--yaml-file YAML_FILE] in batch mode
	 
	 Options:
	    -y, --yaml-file		full path to the yaml file
	    -h, --help			display the help message and exit
	    -v, --verbose		verbose                                                                         
	 
HEREDOC
}

validate()
{
	if [ -f $y ]
	then
		inform " the parameters will be read from "$y
	else
		problem "required yaml file not found! Use the -h option to display the usage"
        exit
	fi
}

function inform
{
    TIMESTAMP=`date +"%a %b %d %X %Y"`
    echo "[$TIMESTAMP] $*" 1>&2
}

function problem
{
    TIMESTAMP=`date +"%a %b %d %X %Y"`
    echo "[$TIMESTAMP] *ERROR*: $*" 1>&2
}

progname=$(basename $0) 

if [[ ! $@ =~ ^\-.+ ]]
then
	echo "No options were passed! See usage below"
	usage
	exit
fi


OPTS=$(getopt -o "hy:v" --long "help,yaml-file:,verbose" -n "$0" -- "$@")
if [ $? != 0 ] ; then echo "Error in command line arguments." >&2 ; usage; exit 1 ; fi
eval set -- "$OPTS"

while true; do
  case "$1" in
    -h | --help ) usage; exit; ;;
    -y | --yaml-file ) y="$2"; shift 2 ;;
    -v | --verbose ) verbose=$((verbose + 1)); shift ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done
validate

workingDirectory=`pwd`


Rscript $workingDirectory/Get_relatedness.R $y $LSB_JOBINDEX $workingDirectory  #$SLURM_ARRAY_TASK_ID    

sleep 1
