#!/bin/bash
find fst -iname '*.fasta' | while read f
do
	d=$(dirname $f)
	test -f $d/raxml.log || \
		bsub "$@" -R "select[tmp>1024],rusage[mem=4096]" -M 4194304 -o $f.out -e $f.err ./tree1k.sh $f 
done
