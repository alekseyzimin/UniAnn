#!/bin/bash
FASTA=$1
HMM_POS=$2
HMM_NEG=$3

function usage {
  echo "Usage: uniann.sh input_fasta positive_hmm_file negative_hmm_file"
}

if [[ -z "$FASTA" ]] || [[ ! -s $FASTA ]]; then
  echo "Fasta file $FASTA not specified or not found!"
  usage
  exit 1
fi
if [[ -z "$HMM_POS" ]] || [[ ! -s $HMM_POS ]]; then
  echo "File with positive HMM models $HMM_POS not specified or not found!"
  usage
  exit 1
fi
if [[ -z "$HMM_NEG" ]] || [[ ! -s $HMM_NEG ]]; then
  echo "File with negative HMM models $HMM_NEG not specified or not found!"
  usage
  exit 1
fi
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
psauron -i $FASTA -a 1>psauron.out 2>&1 && \
$MYPATH/preprocess_scores.pl $FASTA $HMM_POS $HMM_NEG psauron_score.csv &&\
$MYPATH/uniann $FASTA out.ps.txt out.gt.txt out.ag.txt 2>out.err | tee >( grep -v region|gffread -F >$FASTA.gff)
