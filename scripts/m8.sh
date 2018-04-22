#!/bin/bash

#BSUB -L /bin/bash
##BSUB -u iakov.davydov@unil.ch
##BSUB -N
#BSUB -o output-%J-%I.txt
#BSUB -e error-%J-%I.txt
#BSUB -J bs-bur[1-150]
#BSUB –R "select[tmp>1024],rusage[mem=4096]"
#BSUB -M 4194304
##BSUB -n 4
##BSUB -R "span[ptile=4]"

NAME=bur.bin95.filt
SUB=2-m8
SOURCE=genes.95.filt.m0.tgz

source $HOME/mysub/mysub.bash

export OMP_NUM_THREADS=1

GODON=$CLUSTER/godon.06febf8

cmd () {
	if [[ $1 == *.fasta ]]
	then
		fst=$1
		for nwk in ${fst%.*}.*M0.nwk
		do
			out=res/${nwk%.*.*}.M8.log
			tr=${out%.*}.tr
			json=${out%.*}.json
			$GODON test M8 --seed 1 --procs 1 --report 1 --log-level info --no-branch-length --out $out --trajectory $tr --json $json $fst $nwk
		done
	fi
}

run
