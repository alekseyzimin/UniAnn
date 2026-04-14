#!/bin/bash
GC=
RC=
NC=
if tty -s < /dev/fd/1 2> /dev/null; then
    GC='\e[0;32m'
    RC='\e[0;31m'
    NC='\e[0m'
fi

trap abort 1 2 15
function abort {
log "Aborted"
kill -9 0
exit 1
}

log () {
    dddd=$(date)
    echo -e "${GC}[$dddd]${NC} $@"
}

function error_exit {
    dddd=$(date)
    echo -e "${RC}[$dddd]${NC} $1" >&2
    exit "${2:-1}"
}

function usage {
echo "Usage:"
echo "run_gene_finder.sh [arguments]"
echo "-f file sith a single fasta sequence"
echo "-s start PWM models"
echo "-p positive PWM/WAM splice site models"
echo "-n negative PWM/WAM splice site models"

}

#parsing arguments
if [[ $# -eq 0 ]];then
usage
exit 1
fi

while [[ $# > 0 ]]
do
    key="$1"

    case $key in
        -f|--fasta)
            FASTA="$2"
            shift
            ;;
        -s|--start)
            START_PWM="$2"
            shift
            ;;
        -p|--pos)
            POS_PWM="$2";
            shift
            ;;
        -n|--neg)
            NEG_PWM="$2";
            shift
            ;;
        -v|--verbose)
            set -x
            ;;
        -h|--help|-u|--usage)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option $1"
            exit 1        # unknown option
            ;;
    esac
    shift
done

if [[ ! -s $FASTA ]];then
echo "Input file $FASTA not found or not specified!"
usage
exit 1
fi

if [[ ! -s $POS_PWM ]];then
echo "Input file of weights for positive model of splice sites $POS_PWM not found or not specified!"
usage
exit 1
fi

if [[ ! -s $NEG_PWM ]];then
echo "Input file of weights for negative model of splice sites $NEG_PWM not found or not specified!"
usage
exit 1
fi

if [[ ! -s $START_PWM ]];then
echo "Input file of weights for model of start sites $START_PWM not found or not specified!"
usage
exit 1
fi

MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
which psauron
psauron -i $FASTA -a 1>psauron.out 2>&1 && \
$MYPATH/preprocess_psauron_scores.pl $FASTA psauron_score.csv &&\
$MYPATH/compute_markov_scores.pl $FASTA $POS_PWM $NEG_PWM && \
$MYPATH/uniann $FASTA out.ps.txt out.gt.txt out.ag.txt 2>out.err | tee >( grep -v region|gffread -F >$FASTA.gff) > out.txt
echo "Output gff file is $FASTA.gff"
