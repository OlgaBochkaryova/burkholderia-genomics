#!/usr/bin/env python3
from collections import Counter

import numpy as np
from dendropy import Tree

def neg(v):
    assert v in b'01?'
    if v == b'0':
        return b'1'
    elif v == b'1':
        return b'0'
    return b'?'

def negvec(v):
    return np.fromiter((neg(v) for v in v), dtype='S1')

def get_bipart(n, species):
    sub_nodes = set(n.taxon.label for n in n.leaf_iter())
    b1 = np.fromiter((b'1' if l in sub_nodes else b'0'
                      for l in species), dtype='S1')
    b2 = negvec(b1)
    return b1, b2

def compat(bi, b1, b2):
    v1 = np.logical_or(
        bi==b'?', bi==b1
    )
    v2 = np.logical_or(
        bi==b'?', bi==b2
    )
    return np.all(v1) or np.all(v2)

if __name__ == '__main__':
    with open('species') as f:
        species = sorted(l.strip() for l in f)
    tree = Tree.get(path='concat/rooted.nwk', schema='newick')

    bs = []
    bbl = open('bs_branch_lrt.txt', 'w')
    with open('bin70/bs_stat.txt') as f:
        for line in f:
            ds, bi, lrt = line.split()
            lrt = float(lrt)
            bi = np.fromstring(bi, 'S1')
            bs.append((ds, bi, lrt))
    used = set()
    for i, node in enumerate(tree.levelorder_node_iter(lambda n: not n.is_leaf() and n.level() > 0)):
        branch_lrt = []
        nval = 0
        npos = 0
        b1, b2 = get_bipart(node, species)
        print('bootstrap=%s' % node.label)
        try:
            print('edgelength=%0.3g' % node.edge_length, sep='')
        except TypeError:
            print('no edge length, level=%d' % node.level())
        for ds, bi, lrt in bs:
            if compat(bi, b1, b2):
                #assert compat(negvec(bi), b1, b2)
                bi_s = b''.join(bi).decode('utf8')
                v = (ds, bi_s, lrt)
                if v not in used:
                    branch_lrt.append((ds, lrt))
                    used.add(v)
                    nval += 1
                    if lrt > 4:
                        npos += 1
            else:
                #assert not compat(negvec(bi), b1, b2)
                pass
        print('nval=', nval)
        print('npos=', npos)
        print()
        for ds, lrt in branch_lrt:
            print("%d\t%s\t%s" % (i, ds, lrt), file=bbl)
        node.label = '%d' % i

    print(len(used), len(bs))

    ## export all unused
    for ds, bi, lrt in bs:
        bi_s = b''.join(bi).decode('utf8')
        v = (ds, bi_s, lrt)
        if v not in used:
            print("NA\t%s\t%s" % (ds, lrt), file=bbl)
    
    bbl.close()
    
    tree.write(path='bs_tree.nwk', schema='newick')
