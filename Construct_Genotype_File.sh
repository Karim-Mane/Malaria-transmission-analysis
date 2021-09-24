#!/bin/bash

typeset -F SECONDS
module load R
module load bcftools

function usage()
{
	cat << HEREDOC
	
	 About:	This function creates the genotype file from a VCF file. During the processing, the homozygous 
	 reference and alternate alleles will be coded as '0' and '1' respectively. There are 3 ways of recoding the 
	 heterozygous allele: could be considered as a third allele ('raw'), or recoded based into '0' or '1' based
	 on the MAF ('minor_David') or the allelic depth AD ('minor_Karim'). The choice of the method is set using the 
	 -r/--recoding-method option. 
	 
	 Usage:	./$progname [-v INPUT_VCF_FILE] [-o OUTPUT_DIRECTORY] [-m METADATA_FILE] [-r RECODING_METHOD] [-s SELECTIVE_REGIONS_FILE] in interactive mode
	 	bsub -q normal -J "construct_genotype" -n 1 -M5GB -R"select[mem>5GB] rusage[mem=5GB]" -oo construct_genotype.out -eo construct_genotype.err -G team273-vwork 
	 	./$progname [-v INPUT_VCF_FILE] [-o OUTPUT_DIRECTORY] [-m METADATA_FILE] [-r RECODING_METHOD] [-s SELECTIVE_REGIONS_FILE] in batch mode
	 
	 Options:
	    -v, --vcf-file			full path to the input VCF file
	    -o, --output-dir		full path to the output directory
	    -m, --metadata-file		full path to the metadata file
	    -r, --recode-by			heterozygous alleles recoding method
	    -s, --selective-regions			full path to the file with the selective regions; 'none' otherwise
	    -h, --help				display the help message and exit
	    -v, --verbose			verbose                                                                         
	 
HEREDOC
}


function validate()
{
    if [ -f $v ]
    then
        inform "the input VCF file is "$v
    else
        problem " Required input vcf file not found. Use the -h option to display the usage."
        exit
    fi
	
	if [ -d $o ]
    then
        inform $o "the output files will be stored in "$o
    else
        problem " Required output directory not found. Use the -h option to display the usage."
        exit
    fi
	
	if [ -f $m ]
    then
        inform $m "the sample metadata file is "$m
    else
        problem " Required metadata file not found. Use the -h option to display the usage."
        exit
    fi

    if [[ ! -z "$r" ]]   # [[ ! -z "$var" ]] && echo "Not empty"
    then
        inform "the mixed alleles will be recoded using " $r
    else
        problem " the value for the -r option should be of type character. Use the -h option to display the usage."
        exit
    fi
	
	if [ -f $s ] 
    then
        inform $s "the file that contains the selective regions is "$s
	elif [[ ! -z "$s" ]]
	then
		inform "None of the regions of the genome will be removed."
    else
        problem " You must provide a file with selective regions or set the -s option to 'none' if there is no selective regions. 
        Use the -h option to display the usage."
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


OPTS=$(getopt -o "hv:o:m:r:s:v" --long "help,vcf-file:,output-dir:,metadata-file:,recode-by:,selective-regions:,verbose" -n "$0" -- "$@")
if [ $? != 0 ] ; then echo "Error in command line arguments." >&2 ; usage; exit 1 ; fi
eval set -- "$OPTS"

while true; do
  # uncomment the next line to see how shift is working
  # echo "\$1:\"$1\" \$2:\"$2\""
  case "$1" in
    -h | --help ) usage; exit; ;;
    -v | --vcf-file ) v="$2"; shift 2 ;;
    -o | --output-dir ) o="$2"; shift 2 ;;
    -m | --metadata-file ) m="$2"; shift 2 ;;
    -r | --recode-by ) r="$2"; shift 2 ;;
    -s | --selective-regions ) s="$2"; shift 2 ;;
    -v | --verbose ) verbose=$((verbose + 1)); shift ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done

mkdir -p $o
validate

workingDirectory=`pwd`
cd $workingDirectory

Rscript ${workingDirectory}/Construct_Genotype_File.R $v $o $m $s $r

echo -e "\n\nWall time is "$SECONDS" seconds"
