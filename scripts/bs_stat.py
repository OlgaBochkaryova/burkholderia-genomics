#!/usr/bin/env python3
import shelve
import sys
import numpy as np
from dendropy import Tree
from collections import defaultdict


def lrt(j):
    l0 = j['H0']['maxLnL']
    l1 = j['H1']['maxLnL']
    return 2 * l1 - 2 * l0

def neg(v):
    if v == '0':
        return '1'
    elif v == '1':
        return '0'
    elif v == '?':
        return '?'
    else:
        raise ValueError("Unknown state", v)

def get_bipart(ts, species):
    t = Tree.get(data=ts, schema='newick')
    hash_node = t.find_node(lambda n: n.label == '#1')
    sub_nodes = set(n.taxon.label for n in hash_node.leaf_iter())
    all_nodes = set(n.taxon.label for n in t.leaf_node_iter())
    b1 = ''.join('1' if l in sub_nodes
          else ('0' if l in all_nodes else '?')
          for l in species)
    b2 = ''.join(neg(v) for v in b1)
    
    assert b1 != b2
    assert len(b1) == len(b2) and len(b1) == len(species)
    return min(b1, b2)
    

def export_tests(j, d, species):
    if 'tests' in j:
        l = j['tests']
    elif 'tree' in j:
        l = [j] # only a single branch tested, bug in godon
    else:
        return # 0 branches tested
    for t in l:
        bp = get_bipart(t['tree'], species)
        d[bp].append((lrt(t), t['H1']['maxLParameters']['omega2']))
    

if __name__ == '__main__':
    with open('species') as f:
        species = sorted(l.strip() for l in f)
    db = shelve.open('bur.db')
    files = [k for k in db.keys() if k != 'datasets']
    files.sort()
    ds2f = defaultdict(list)
    for fn in files:
        try:
            dataset, bs, model, _ = fn.split('.')
        except ValueError:
            dataset, model, _ = fn.split('.')
            bs = 'ml'
        if model == 'BS':
            ds2f[dataset].append(fn)
        
    for i, ds in enumerate(ds2f):
        if i % 10 == 0:
            print('%d/%d (%.0f%%)' % (i, len(ds2f), i/len(ds2f)*100), file=sys.stderr)
        tests = defaultdict(list)
        for fn in ds2f[ds]:
            export_tests(db[fn], tests, species)
        for bp, l in tests.items():
            if len(l) > 3:
                ml = min(l)
                print(ds, bp, ml[0], ml[1])
