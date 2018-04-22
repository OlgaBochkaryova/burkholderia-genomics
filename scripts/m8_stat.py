#!/usr/bin/env python3
import shelve
import numpy as np

def lrt(j):
    l0 = j['H0']['maxLnL']
    l1 = j['H1']['maxLnL']
    return 2 * l1 - 2 * l0
    

if __name__ == '__main__':
    db = shelve.open('bur.db')
    files = [k for k in db.keys() if k != 'datasets']
    files.sort()
    for ds in db['datasets']:
        l = []
        for fn in files:
            try:
                dataset, bs, model, _ = fn.split('.')
            except ValueError:
                dataset, model, _ = fn.split('.')
                bs = 'ml'
            if dataset == ds and model == 'M8':
                l.append((lrt(db[fn]), db[fn]['H1']['maxLParameters']['omega']))
        ml = min(l)
        print(ds, ml[0], ml[1])
