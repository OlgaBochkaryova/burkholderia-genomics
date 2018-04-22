#!/usr/bin/env python3
import json
import os
import shelve
import glob
import sys
from collections import defaultdict

if __name__ == '__main__':
    datasets = set()
    db = shelve.open('bur.db')
    files = glob.glob('res/*.json')
    for i, fn in enumerate(files):
        bfn = os.path.basename(fn)
        ds = bfn.split('.')[0]
        datasets.add(ds)
        if i%100 == 0:
            print('%d/%d (%0.0f%%)' % (i, len(files), float(i)/len(files)*100), file=sys.stdout, end='\r')
            sys.stdout.flush()
        with open(fn) as f:
            j = json.load(f)
            db[bfn] = j
    print('done')
    db['datasets'] = datasets
    db.close()
