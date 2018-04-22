#!/usr/bin/env python3
import glob
import os.path

from dendropy import Tree

print('OG', 'tlen')
for fn in glob.iglob('binning.fulltree/*/RAxML_bestTree.tree'):
    t = Tree.get(path=fn, schema='newick')
    b = os.path.basename(os.path.dirname(fn))
    print(b, t.length())
