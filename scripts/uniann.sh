#!/bin/bash
FASTA=$1
HMM_POS=$2
HMM_NEG=$3
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
psauron -i $FASTA -a 1>psauron.out 2>&1 && \
$MYPATH/preprocess_scores.pl $FASTA $HMM_POS $HMM_NEG psauron_score.csv &&\
$MYPATH/uniann $FASTA out.ps.txt out.gt.txt out.ag.txt 2>out.err | tee >( grep -v region|gffread -F >$FASTA.gff)
