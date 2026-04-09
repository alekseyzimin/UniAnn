#!/bin/bash
VERSION=1.0.0
echo "Making distribution" && \
rm -rf uniann && \
mkdir -p uniann && \
cp -r LICENSE install.sh gffcompare gffread src scripts uniann && \
tar cvzf uniann.$VERSION.tgz uniann && \
rm -rf uniann && \
echo "Done"
