#!/bin/bash
set -xe
ROOT=$PWD && \
mkdir -p $ROOT/bin && \
(cd src && make && cp uniann $ROOT/bin) && \
(cd scripts && cp *.pl *.sh $ROOT/bin && chmod 0755 $ROOT/bin/*.{pl,sh}) && \
echo "All done" || echo "Installation failed"
