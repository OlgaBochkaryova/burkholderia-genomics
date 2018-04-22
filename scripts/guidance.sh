#!/bin/bash
fn=$(readlink -f $1)
PATH=/home/idavydov/mafft/mafft-linux64:$PATH
outdir=$fn.out
mkdir -p $outdir
perl ~/software/guidance/guidance.v2.01/www/Guidance/guidance.pl --seqFile $fn --msaProgram MAFFT --seqType codon --outDir $outdir
