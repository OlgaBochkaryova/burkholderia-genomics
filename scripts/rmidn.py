#!/usr/bin/env python
import glob
from collections import defaultdict
from Bio import SeqIO, AlignIO

## this function not used; we keep gaps in the ali
def rmgaps(ali):
    bad_cols = set()
    for rec in ali:
        for i, nuc in enumerate(rec.seq):
            if nuc.upper() in ('X', 'N', '-'):
                bad_cols.add(i)
    print sorted(bad_cols)
    
    ranges = []
    start = None
    last = None
    for col in sorted(bad_cols):
        if start is None:
            start = col
        else:
            if col != last + 1:
               ranges.append((start, last + 1))
               start = col
        last = col
    if last is not None:
        ranges.append((start, last + 1))

    res = ali[:, :0]
    start = 0
    for st, en in ranges:
        res += ali[:, start:st]
        start = en
    res += ali[:, start:ali.get_alignment_length()]
    assert res.get_alignment_length() == ali.get_alignment_length() - len(bad_cols)
    return res

for fn in glob.iglob('*/.fasta'):
    ali = AlignIO.read(fn, 'fasta')
    used_seq = set()
    new_recs = []
    seq2rec = defaultdict(list)
    for rec in ali:
        seq = str(rec.seq)
        seq2rec[seq].append(rec)
    for recs in seq2rec.values():
        new_recs.append(min(recs, key=lambda r: r.id))
    if len(ali) != len(new_recs):
        print 'removed %d seqs' % (len(ali) - len(new_recs))
    nfn = fn.rsplit('.', 1)[0] + '.flt.fasta'
    SeqIO.write(new_recs, nfn, 'fasta')
