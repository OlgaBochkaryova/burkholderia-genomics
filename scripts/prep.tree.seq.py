#!/usr/bin/env python3
## this script takse supergenes, and
## prepares data for positive selection analysis
import sys
import os
import os.path
import glob
import random
import numpy as np

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from dendropy import Tree, TreeList

def remove_identical(t, ali):  
    used_seqs = set()
    acs_to_rm = set()
    seqs = {record.id: str(record.seq) for record in ali}

    for node in t.leaf_node_iter():
        ac = node.taxon.label
        if ac not in seqs or seqs[ac] in used_seqs:
            acs_to_rm.add(ac)
        else:
            used_seqs.add(seqs[ac])
    orig_leafs = [node.taxon.label for node in t.leaf_node_iter()]
    t.prune_taxa_with_labels(acs_to_rm)
    upd_nleafs = [node.taxon.label for node in t.leaf_node_iter()]
    assert len(orig_leafs) - len(upd_nleafs) == len(acs_to_rm)
    return acs_to_rm

def remove_extra_seq(t, ali):
    in_ali = {record.id for record in ali}
    in_tree = {tip.taxon.label for tip in t.leaf_node_iter()}
    bad = in_ali - in_tree
    return remove_bad(ali, bad)

def remove_bad(ali, bad):
    return MultipleSeqAlignment([record for record in ali if not
                                 record.id in bad])

def tree_upper(t):
    for taxon in t.taxon_namespace:
        taxon.label = taxon.label.upper()
    return t

def ali_upper(ali):
    for record in ali:
        record.id = record.id.upper()
        record.description = ''
    return ali

def tree_topo(t):
    return t.as_string(
        suppress_internal_node_labels=True,
        unquoted_underscores=True,
        suppress_edge_lengths=True,
        schema='newick')

def equal_branch_lengths(t, blen):
    for node in t.postorder_node_iter():
        if node.level() == 0: #root
            continue
        node.edge_length = blen

class Bins:
    def __init__(self, path):
        self.path = path

    def read_bins(self):
        self.bin2genes = {}
        for fn in sorted(glob.iglob(os.path.join(self.path, 'pair', 'bin.*.txt'))):
            self.bin2genes[os.path.basename(fn)] = [l.strip() for l in open(fn)]

    def read_gene(self, gene):
        return ali_upper(AlignIO.read(os.path.join(self.path, 'binning.fulltree', gene, gene + '.fasta'),
                            'fasta'))

    def read_genes_for_bin(self, bin_name):
        return {gene: self.read_gene(gene) for gene in self.bin2genes[bin_name]}

    def get_ml_tree(self, bin_name):
        return tree_upper(Tree.get(
            path=os.path.join(
                self.path, 'supergenes',
                bin_name, 'RAxML_bipartitions.bipart'),
            preserve_underscores=True,
            schema='newick'))

    def get_bs_trees(self, bin_name):
        tl = TreeList.get(
            path=os.path.join(
                self.path, 'supergenes',
                bin_name, 'RAxML_bootstrap.bootstrap'),
            preserve_underscores=True,
            schema='newick')
        tree_upper(tl[0])
        return tl

    def save_atree(self, t, a, gene):
        bfn = os.path.join(self.path, 'genes', gene)
        t.write(path=bfn + '.nwk', 
                   unquoted_underscores=True,
                   schema='newick')
        f = open(bfn + '.fasta', 'w')
        AlignIO.write(a, f, 'fasta')
        f.close()

    def save_tree(self, t, gene, i):
        tfn = os.path.join(self.path, 'genes', gene + '.%03d.nwkb' % i)
        t.write(path=tfn, 
                   unquoted_underscores=True,
                   schema='newick')

    def process_bin(self, bin_name):
        genes = self.read_genes_for_bin(bin_name)
        mltree = self.get_ml_tree(bin_name)
        mblen = np.mean([node.edge_length for node in mltree.preorder_node_iter() if node.edge_length and node.level() > 0])
        bstrees = self.get_bs_trees(bin_name)
        print(bin_name)
        for gene, ali in genes.items():
            t = mltree.clone(depth=1)
            removed_leafs = remove_identical(t, ali)
            a = remove_extra_seq(t, ali)
            assert len(a) == len(list(t.leaf_node_iter()))
            self.save_atree(t, a, gene)

            ## now save all the bootstrap trees
            trees = []
            for bt in bstrees[:5]:
                t = bt.clone(depth=1)
                t.prune_taxa_with_labels(removed_leafs)
                trees.append(t)
            ## export trees
            for i, t in enumerate(trees):
                equal_branch_lengths(t, mblen)
                self.save_tree(t, gene, i)
                assert len(list(t.leaf_node_iter())) == len(a)
                    

    def process(self):
        for bin_name in self.bin2genes:
            self.process_bin(bin_name)
            
        
if __name__ == '__main__':
    random.seed(1)
    path = sys.argv[1]
    bins = Bins(path)
    bins.read_bins()
    bins.process()
