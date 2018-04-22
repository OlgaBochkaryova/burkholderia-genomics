#!/bin/bash
d=$(dirname $(readlink -f $1))
f=$(basename $1)
cd $d
log=raxml.log
test -f $log && exit 0
rm -f $log
raxml=${RAXML:-$HOME/raxml/standard-RAxML/raxmlHPC-PTHREADS-AVX}
$raxml -m GTRGAMMA -T 10 -s $f -n tree -N 20 -p 1 2>&1 >> $log && \
$raxml -m GTRGAMMA -T 10 -s $f -n bootstrap -N 100 -p 1 -b 1 2>&1 >> $log && \
$raxml -m GTRGAMMA -f b -t RAxML_bestTree.tree -z RAxML_bootstrap.bootstrap -n bipart 2>&1 >> $log

if [ ! $? -eq 0 ]
then
    echo error in $1
fi
exit 0
