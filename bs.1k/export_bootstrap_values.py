#!/usr/bin/env python3
import sys
from dendropy import Tree
from glob import glob
from os.path import dirname, basename, join

def export_tree(n, fn):
    t = Tree.get(path=fn, schema='newick')

def tostr(t):
    return t.as_string(schema='newick', suppress_internal_node_labels=True)

def get_node_labels(t):
    for node in t.internal_nodes():
        if node.label:
            yield node.label

def main():
    d1 = sys.argv[1]
    d2 = sys.argv[2]
    d1_name = basename(d1)
    d2_name = basename(d2)
    print('og {} {}'.format(d1_name, d2_name))
    d1_files = list(sorted(glob(join(d1, '*', 'RAxML_bipartitions.bipart'))))
    d2_files = list(sorted(glob(join(d2, '*', 'RAxML_bipartitions.bipart'))))
    assert len(d1_files) == len(d2_files)
    for fn1, fn2 in zip(d1_files, d2_files):
        t1 = Tree.get(path=fn1, schema='newick')
        t2 = Tree.get(path=fn2, schema='newick')
        assert tostr(t1) == tostr(t2)
        t1_og = basename(dirname(fn1))
        t2_og = basename(dirname(fn2))
        assert t1_og == t2_og
        labs = zip(get_node_labels(t1), get_node_labels(t2))
        for l1, l2 in labs:
            print(t1_og, l1, l2)

if __name__ == '__main__':
    main()
