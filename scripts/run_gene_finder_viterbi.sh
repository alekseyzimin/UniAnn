#!/bin/bash
FASTA=$1
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
echo "This code only works for Drosophila"
psauron -i $FASTA -a 1>psauron.out 2>&1 && \
$MYPATH/preprocess_scores.pl $FASTA $MYPATH/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.pwm $MYPATH/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.neg.pwm psauron_score.csv &&\
$MYPATH/gene_finder_viterbi.pl $FASTA out.ps.txt out.gt.txt out.ag.txt 2>out.err | tee >( grep -v region|gffread -F >$FASTA.gff)
