#!/bin/bash
set -o pipefail
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"

FASTA="genome.fa"
CODING_FRAME="coding.txt"
SCOREFILE="scores.txt"
MULT=`perl -e 'print exp(1)'`
OUTDEV="out.err"
MIN_CDS=200
BOTH_STRANDS=0
EXTRA_FLAGS=()

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
echo "-s file with scores for starts, stops and splice sites, 1-based coordinates"
echo "-c coding frame score file"
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
        -c|--coding)
            CODING_FRAME="$2"
            shift
            ;;
        -m|--mult)
            MULT="$2"
            shift
            ;;
        -n|--noviterbi)
            OUTDEV="/dev/null"
            EXTRA_FLAGS+=(-n)
            ;;
        -s|--scores)
            SCOREFILE="$2"
            shift
            ;;
        -v|--verbose)
            set -x
            EXTRA_FLAGS+=(-v)
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
echo "Input training annotation file $SCOREFILE not found or not specified!"
usage
exit 1
fi

if [[ ! -s $CODING_FRAME ]];then
echo "Input file of coding frame scores $CODING_FRAME not found or not specified!"
usage
exit 1
fi

# For -a, the helper prepares separate inputs and runs this forward pipeline for each strand.
# It maps reverse predictions back to the original FASTA coordinates and combines both GFFs.
# this need rework
#if [[ $BOTH_STRANDS == 1 ]]; then
#    exec perl "$MYPATH/uniann_both_strands.pl" "$MYPATH/uniann.sh" \
#        "$FASTA" "$PSAURON" "$SCOREFILE" "$MULT" "${EXTRA_FLAGS[@]}"
#fi

#this produces out.ps.txt
log "Preprocessing inputs" && \
"$MYPATH/convert_coding_frame.pl" "$FASTA" < "$CODING_FRAME" > out.ps.txt.tmp && \
mv out.ps.txt.tmp out.ps.txt || error_exit "Failed to prepare coding frame scores" "$?"

FACTOR=$(perl -ane '{$F[5]=$F[6] if($#F>5);;print join("\t",@F),"\n" if($F[2] eq "+");}' "$SCOREFILE" | \
  perl -ane 'BEGIN{$max_score=0}{if($F[3] eq "donor"){$score=log($F[5]*'"$MULT"'+1e-10);$max_score=$score if($score>$max_score);}}END{die("Incorrect scores in the input file: must be between 0 and 1!") if($max_score<=0);print int(1000/$max_score+0.5)}') || \
  error_exit "Failed to calculate score factor" "$?"
log "Multiplier is $MULT, Factor is $FACTOR"

# Write each site file synchronously so failures stop before decoding.
for SITE in donor acceptor start stop; do
    case "$SITE" in
        donor) SITE_FILE=out.gt.txt ;;
        acceptor) SITE_FILE=out.ag.txt ;;
        start) SITE_FILE=out.atg.txt ;;
        stop) SITE_FILE=out.stop.txt ;;
    esac
    perl -ane '{$F[5]=$F[-1] if($#F>5);;print join("\t",@F),"\n" if($F[2] eq "+");}' "$SCOREFILE" | \
      perl -ane '{if($F[3] eq "'"$SITE"'"){$score=log($F[5]*('"$MULT"')+1e-10)*'"$FACTOR"'; print $F[1]-1,"\t$score\n"}}END{close STDOUT or die "Cannot write site scores: $!"}' > "$SITE_FILE" || \
      error_exit "Failed to prepare $SITE scores" "$?"
done

log "Building gene models" && \
#also enforce MIN_CDS
$MYPATH/uniann $FASTA out.ps.txt out.gt.txt out.ag.txt out.atg.txt out.stop.txt 2>$OUTDEV | \
  gffread --tlf |\
  perl -F'\t' -ane '{
    if($F[8]=~/exonCount=(1|2);exons=(\S+);CDS=(\d+):(\d+);CDSphase=\d/){
      $cdslen=0;
      @exons=split(/,/,$2);
      foreach $e(@exons){
        $cdslen+=$2-$1+1 if($e=~/(\d+)-(\d+)/);
      }
      print if($cdslen>'$MIN_CDS');
    }else{
      print;
    }
  }' |\
  gffread > "$FASTA.gff.tmp" || error_exit "Failed to build gene models" "$?"
mv "$FASTA.gff.tmp" "$FASTA.uniann.gff" || error_exit "Failed to save output GFF" "$?"

echo "Output gff file is $FASTA.uniann.gff"
