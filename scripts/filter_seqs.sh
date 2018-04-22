#!/bin/bash
P=/home/idavydov/software/guidance/guidance.v2.01/www/Guidance/
thr=0.93
i=0
N=1
for fn in ogs/*.fas
do
    base=$fn.out
    if [ -f $base/MSA.MAFFT.aln.With_Names ]
    then
	((i++))
	in=$base/MSA.MAFFT.aln.With_Names
	out=$in.ResFilt
	echo perl $P/maskLowScoreResidues.pl $in $base/MSA.MAFFT.Guidance2_res_pair_res.scr $out $thr nuc
	perl $P/maskLowScoreResidues.pl $in $base/MSA.MAFFT.Guidance2_res_pair_res.scr $out $thr nuc
	in=$out
	out=$in.SeqFilt
	echo perl $P/Remove_Seq_bellow_Cutoff.pl --MSA $in --Scores $base/MSA.MAFFT.Guidance2_res_pair_seq.scr --FilterdSeq $out --Cutoff 0.8 --RemovedSeq $out.Removed --Type ByRowNum 
	perl $P/Remove_Seq_bellow_Cutoff.pl --MSA $in --Scores $base/MSA.MAFFT.Guidance2_res_pair_seq.scr --FilterdSeq $out --Cutoff 0.8 --RemovedSeq $out.Removed --Type ByRowNum 
	echo $i $fn
    fi
done
