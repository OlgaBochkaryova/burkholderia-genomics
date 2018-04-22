#!/usr/bin/env python3
import glob
import os
import numpy
from Bio import SeqIO

def avggc(seqs):
    at = len([l for rec in seqs for l in rec.seq if l in 'AT'])
    gc = len([l for rec in seqs for l in rec.seq if l in 'GC'])
    if at > 0 and gc > 0:
        return float(gc)/(gc+at)

def deltagc(seqs):
    gcs = []
    for rec in seqs:
        at = len([l for l in rec.seq if l in 'AT'])
        gc = len([l for l in rec.seq if l in 'GC'])
        if at > 0 and gc > 0:
            gcs.append(float(gc)/(gc+at))
    return numpy.std(gcs)

print('OG', 'GC', 'GC.sd')
for fn in glob.iglob('genes/*.fasta'):
    seqs = list(SeqIO.parse(fn, 'fasta'))
    print(os.path.basename(fn).rsplit('.')[0], avggc(seqs), deltagc(seqs))
