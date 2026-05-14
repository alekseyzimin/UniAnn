#!/bin/bash

FASTA="genome.fa"
TRAINING_G="training_genome.fa"
TRAINING_A="training_annotation.gff"
PSAURON="psauron_score.csv"
MIN_CDS=300
MIN_SINGLE_CDS=400

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
echo "MUST HAVE out.gt.txt and out.ag.txt splice site score files in the folder!!!"
echo "uniann.sh [arguments]"
echo "-f file sith a single fasta sequence"
echo "-g genome to train on"
echo "-a annotation to train on"
echo "-p psauron score file"
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
        -a|--annotation)
            TRAINING_A="$2"
            shift
            ;;
        -g|--genome)
            TRAINING_G="$2"
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

if [[ ! -s $TRAINING_G ]];then
echo "Input training genome file $TRAINING_G not found or not specified!"
usage
exit 1
fi

if [[ ! -s $TRAINING_A ]];then
echo "Input training annotation file $TRAINING_A not found or not specified!"
usage
exit 1
fi

if [[ ! -s $PSAURON ]];then
echo "Input file of PSAURON scores $PSAURON not found or not specified!"
usage
exit 1
fi

MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"

#this produces out.ps.txt
log "Preprocessing psauron scores" && \
$MYPATH/preprocess_psauron_scores.pl $FASTA $PSAURON && \
#this produces out.atg.txt
log "Scoring candidate start sites" && \
$MYPATH/compute_start_scores.pl $TRAINING_G $TRAINING_A > start.pwm && \
$MYPATH/compute_stop_scores.pl $TRAINING_G $TRAINING_A > stop.pwm && \
$MYPATH/score_start_sites.pl start.pwm < $FASTA && \
$MYPATH/score_stop_sites.pl stop.pwm < $FASTA && \
#this produces out.gt.txt and out.ag.txt
#log "Scoring candidate splice sites" && \
#$MYPATH/compute_markov_scores $FASTA $POS_PWM $NEG_PWM && \
log "Building gene models" && \
$MYPATH/uniann $FASTA out.ps.txt out.gt.txt out.ag.txt out.atg.txt out.stop.txt 2>out.err| \
  tee >( perl -ane '{
          push @lines,$_;
        }END{
          $prev="region";
          for(my $i=0;$i<$#lines;$i++){
            @f=split(/\t/,$lines[$i]);
            @ff=split(/\t/,$lines[$i+1]);
            if($prev eq "region" && $f[2] eq "CDS" && $ff[2] eq "region"){
            #AVERAGE SCORE MUST BE ABOVE .75 FOR SINGLE EXON
              print $lines[$i] if($f[5]>0.75);
            }else{
              print $lines[$i] unless($f[2] =~ /region|intron/);
            }
            $prev=$f[2]; 
          }
          @f=split(/\t/,$lines[-1]);
          print $lines[-1] unless($f[2] =~ /region|intron/);
        }'  |\
    gffread --tlf | \
    perl -F'\t' -ane '{
      if($F[8]=~/ID=(\S+);exonCount=(\d+);exons=(\S+);CDS=(\d+):(\d+);CDSphase/){
        print if(($2==1 && $5-$4 > '$MIN_SINGLE_CDS') || ($2 > 1 && $5-$4 > '$MIN_CDS'));
      }
    }' |\
    gffread -F  >$FASTA.gff.tmp) > $FASTA.all.gff.tmp && \
mv $FASTA.gff.tmp $FASTA.gff && \
mv $FASTA.all.gff.tmp $FASTA.all.gff && \
echo "Output gff file is $FASTA.gff"
