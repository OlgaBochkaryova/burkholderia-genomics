#!/bin/bash
d=$(dirname $(readlink -f $1))
f=$(basename $1)
raxml=$(readlink -f ./raxmlHPC-SSE3)
cd $d
log=raxml.log
test -f $log && exit 0
rm -f $log
$raxml -m GTRGAMMA -s $f -n bootstrap -N 1000 -p 1 -b 1 2>&1 >> $log

if [ ! $? -eq 0 ]
then
    echo error in $1
fi
exit 0
