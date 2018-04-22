#!/usr/bin/env python3
import sys
import glob
import os.path
from collections import defaultdict

from Bio.SeqIO.FastaIO import SimpleFastaParser

def wrap(st, n=70):
    return (st[i:i+n] for i in range(0, len(st), n))

if __name__ == '__main__':
    with open('species') as f:
        species = set(l.strip() for l in f)
    seqs = defaultdict(list)
    for path in glob.iglob('ogs/*.out'):
        if os.path.getsize(os.path.join(path, 'MSA.MAFFT.aln.With_Names.ResFilt.SeqFilt.Removed')) > 0:
            print('bad sequence', path, file=sys.stderr)
            continue
        with open(os.path.join(path, 'MSA.MAFFT.aln.With_Names')) as f:
            for i, (rid, seq) in enumerate(SimpleFastaParser(f)):
                  seqs[rid].append(seq)
        assert i == len(species) - 1

    for rid, seql in seqs.items():
        print('>', rid, sep='')
        for l in wrap(''.join(seql)):
            print(l)
        
