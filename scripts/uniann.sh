#!/bin/bash
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"

FASTA="genome.fa"
PSAURON="psauron_score.csv"
SCOREFILE="scores.txt"
MULT=`perl -e 'print exp(1)'`
VITERBI=1

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
echo "uniann.sh [arguments]"
echo "-f file sith a single fasta sequence"
echo "-m multiplier to rescale probabilities before applying log, must be a number between 1 and 100, default: exp(1)"
echo "-s file with AI-derived scores"
echo "-p psauron score file"
echo "-n (flag) do not save Viterbi matrix in out.err"
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
        -p|--psauron)
            PSAURON="$2"
            shift
            ;;
        -m|--mult)
            MULT="$2"
            shift
            ;;
        -n|--noviterbi)
            VITERBI=0
            shift
            ;;
        -s|--scores)
            SCOREFILE="$2"
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

if [[ ! -s $SCOREFILE ]];then
echo "Input training annotation file $TRAINING_A not found or not specified!"
usage
exit 1
fi

if [[ ! -s $PSAURON ]];then
echo "Input file of PSAURON scores $PSAURON not found or not specified!"
usage
exit 1
fi

#this produces out.ps.txt
log "Preprocessing inputs" && \
if [ ! -s out.ps.txt ];then
  $MYPATH/preprocess_psauron_scores.pl $FASTA $PSAURON
fi

FACTOR=`cat $SCOREFILE |perl -ane '{$F[5]=$F[6] if($#F>5);;print join("\t",@F),"\n" if($F[2] eq "+");}'|perl -ane 'BEGIN{$max_score=0}{if($F[3] eq "donor"){$score=log($F[5]*'$MULT'+1e-10);$max_score=$score if($score>$max_score);}}END{die("Incorrect scores in the input file: must be between 0 and 1!") if($max_score<=0);print int(1000/$max_score+0.5)}'` && \
log "Multiplier is $MULT, FACTOR is $FACTOR" && \
#this produces out.atg.txt out.gt.txt and out.ag.txt out.stop
cat $SCOREFILE | \
  perl -ane '{$F[5]=$F[-1] if($#F>5);;print join("\t",@F),"\n" if($F[2] eq "+");}' | \
  tee >(perl -ane '{if($F[3] eq "donor"){$score=log($F[5]*'$MULT'+1e-10)*'$FACTOR'; print $F[1]-1,"\t$score\n"}}' > out.gt.txt) |\
  tee >(perl -ane '{if($F[3] eq "acceptor"){$score=log($F[5]*'$MULT'+1e-10)*'$FACTOR'; print $F[1]-1,"\t$score\n"}}' > out.ag.txt) |\
  tee >(perl -ane '{if($F[3] eq "start"){$score=log($F[5]*('$MULT')+1e-10)*'$FACTOR'; print $F[1]-1,"\t$score\n"}}' > out.atg.txt) | \
  perl -ane '{if($F[3] eq "stop"){$score=log($F[5]*('$MULT')+1e-10)*'$FACTOR'; print $F[1]-1,"\t$score\n"}}' > out.stop.txt && \

log "Building gene models" && \
if [ $VITERBI -lt 1 ];then
  $MYPATH/uniann $FASTA out.ps.txt out.gt.txt out.ag.txt out.atg.txt out.stop.txt 2>/dev/null|gffread -F > $FASTA.gff.tmp && \
  mv $FASTA.gff.tmp $FASTA.uniann.gff
else
  $MYPATH/uniann $FASTA out.ps.txt out.gt.txt out.ag.txt out.atg.txt out.stop.txt 2>out.err|gffread -F > $FASTA.gff.tmp && \
  mv $FASTA.gff.tmp $FASTA.uniann.gff
fi 

echo "Output gff file is $FASTA.uniann.gff"
