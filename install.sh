#!/bin/bash
set -xe
ROOT=$PWD && \
mkdir -p $ROOT/bin && \
cp gffcompare gffread $ROOT/bin && \
(cd src && make && mv uniann compute_markov_scores $ROOT/bin) && \
(cd scripts && cp *.pl *.sh $ROOT/bin && chmod 0755 $ROOT/bin/*.{pl,sh}) && \
echo "All done" || echo "Installation failed"
