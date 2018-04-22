#!/usr/bin/env python3
import shelve
import numpy as np

def get_par(j):
    res = {}
    l0 = j['H0']['maxLnL']
    l1 = j['H1']['maxLnL']
    res['lrt'] = 2 * l1 - 2 * l0
    res.update(j['H1']['maxLParameters'])
    return res

if __name__ == '__main__':
    print('OG', 'tree', 'LRT', 'p', 'q', 'omega', 'kappa', 'p0')
    db = shelve.open('bur.db')
    for fn in db:
        if fn.endswith('.json'):
            try:
                dataset, bs, model, _ = fn.split('.')
            except ValueError:
                dataset, model, _ = fn.split('.')
                bs = 'ml'
            if model == 'M8':
                par = get_par(db[fn])
                print(dataset, bs, par['lrt'],
                      par['p'], 
                      par['q'], 
                      par['omega'], 
                      par['kappa'], 
                      par['p0']
                )
