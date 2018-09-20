#!/bin/bash
d="$(dirname $1)"
cd $d
log=raxml.log
test -f $log && exit 0
rm -f $log
raxml=${RAXML:-$HOME/raxml/standard-RAxML/raxmlHPC-PTHREADS-AVX}
$raxml -m GTRGAMMA -f b -t RAxML_bestTree.tree -z RAxML_bootstrap.bootstrap -n bipart 2>&1 >> $log

if [ ! $? -eq 0 ]
then
    echo error in $1
fi
exit 0
