#!/bin/bash

#BSUB -L /bin/bash
##BSUB -u iakov.davydov@unil.ch
##BSUB -N
#BSUB -o output-%J-%I.txt
#BSUB -e error-%J-%I.txt
#BSUB -J bs-bur[1-80]
#BSUB –R "select[tmp>1024],rusage[mem=4096]"
#BSUB -M 4194304
##BSUB -n 4
##BSUB -R "span[ptile=4]"

NAME=bur.bin95.filt
SUB=1-bs
SOURCE=genes.95.filt.tgz

source $HOME/mysub/mysub.bash

export OMP_NUM_THREADS=1

GODON=$CLUSTER/godon.06febf8

cmd () {
	if [[ $1 == *.fasta ]]
	then
		fst=$1
		nwk=${1%.*}.nwk
		out=res/${fst%.*}.BS.log
		tr=${out%.*}.tr
		json=${out%.*}.json
		$GODON test --seed 1 --procs 1 --no-leaves --log-level info --m0-tree --out $out --trajectory $tr --json $json BS $fst $nwk
		for nwkb in ${fst%.*}.*.nwkb
		do
			nwk=$nwkb
			out=res/${nwk%.*}.BS.log
			tr=${out%.*}.tr
			json=${out%.*}.json
			$GODON test --seed 1 --procs 1 --no-leaves --log-level info --m0-tree --out $out --trajectory $tr --json $json BS $fst $nwk
		done
	fi
}

run
