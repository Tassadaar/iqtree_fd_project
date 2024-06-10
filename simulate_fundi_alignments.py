#!/usr/bin/env python

# standard modules
import argparse
import re
import sys
import subprocess

# stuff you need to install
from ete3 import Tree
from Bio import AlignIO

# argparse
parser = argparse.ArgumentParser(
        description=
        '''
        Functional Divergence Alignment Simulator
        ''')
parser.add_argument("-te", "--tree", 
                    required=True, 
                    help="Tree file in newick format, must be rooted")
parser.add_argument("-d", "--definition", 
                    required=True,
                    help="Definition file that splits the tree by FunDi branch")
parser.add_argument('-m', '--model',
                   required=False,
                   default='LG',
                   help='Model to simulate alignment under')
parser.add_argument('-l', '--length',
                    required=True,
                    type=int,
                    help='Total simulated alignment length')
parser.add_argument('-r', '--rho',
                    required=True,
                    type=float,
                    help='Rho: proportion of sites to have Functional Divergence in the simulated alignment')
parser.add_argument('-n', '--replicates',
                    required=False,
                    default=5,
                    type=int,
                    help='Number of replicate simulations to generate')

args = parser.parse_args()




def main(args):

    # load full tree
    master_tree = Tree(args.tree)
    # print(master_tree)

    # load and validate definition file
    defined_groups = load_validate_def_file(args.definition, master_tree)
    # print(defined_groups)

    # split tree based on definition file
    ab_tree, a_tree, b_tree = get_simulation_trees(master_tree, defined_groups)
    # print(ab_tree, a_tree, b_tree)
    # write subtrees into newick files
    ab_tree.write(format=1, outfile='T_ab.newick')
    a_tree.write(format=1, outfile="T_a.newick")
    b_tree.write(format=1, outfile="T_b.newick")

    for n in range(1, args.replicates+1):

        # simulate (1-rho) fraction of alignment with full tree
        print(f'Running AliSim for T_ab alignment, replicate {n:02d}', file=sys.stderr)
        alisim_command = [
                'iqtree2',
                '--alisim', f'A_ab_{n:02d}',
                '-t', 'T_ab.newick',
                '-m', args.model,
                '--length', str( round((1-args.rho)*args.length) ),
                '-a', '0.5',
                '--seqtype', 'AA',
                '-af', 'fasta',
                '-nt', '4',
                '--quiet'
                ]
        # print(alisim_command)
        subprocess.run(alisim_command)

        # simulate (rho) fraction of alignment with subtree 1
        print(f'Running AliSim for T_a alignment, replicate {n:02d}', file=sys.stderr)
        alisim_command = [
                'iqtree2',
                '--alisim', f'A_a_{n:02d}',
                '-t', 'T_a.newick',
                '-m', args.model,
                '--length', str( round(args.rho*args.length) ),
                '-a', '0.5',
                '--seqtype', 'AA',
                '-af', 'fasta',
                '-nt', '4',
                '--quiet'
                ]
        subprocess.run(alisim_command)

        # simulate (rho) fraction of alignment with subtree 2
        print(f'Running AliSim for T_b alignment, replicate {n:02d}', file=sys.stderr)
        alisim_command = [
                'iqtree2',
                '--alisim', f'A_b_{n:02d}',
                '-t', 'T_b.newick',
                '-m', args.model,
                '--length', str( round(args.rho*args.length) ),
                '-a', '0.5',
                '--seqtype', 'AA',
                '-af', 'fasta',
                '-nt', '4',
                '--quiet'
                ]
        subprocess.run(alisim_command)

        # merge 3 sub alignments into single full alignment

        # first load simulated sub_alignments
        A_ab = AlignIO.read(f'A_ab_{n:02d}.fa', 'fasta')
        A_a  = AlignIO.read(f'A_a_{n:02d}.fa', 'fasta')
        A_b  = AlignIO.read(f'A_b_{n:02d}.fa', 'fasta')

        # initiate an empty MultipleSeqAlignment object to store 
        # merged alignment in
        A_a_A_b = AlignIO.MultipleSeqAlignment([])

        # use .extend() to merge vertically A_a and A_b
        A_a_A_b.extend(A_a)
        A_a_A_b.extend(A_b)

        # finally merge left and right
        # (ensure that sequence order is consistent between the two alignments!)
        A = A_ab + A_a_A_b
        AlignIO.write(A, f'A_{n:02d}.aln', 'fasta')

def load_validate_def_file(def_file, tree):
    # store definition file in memory
    # as a list of two sets of taxa
    # this ensures uniqueness, assuming that order does not matter with def files
    with open(def_file, "r") as file:
        taxa_groups = []
        for line in file:
            # clean line by removing whitespace, comma or newline characters
            clean_line = line.rstrip(" ,\n")
            # ignore empty lines
            if clean_line == "": continue
            # parse clean line delimited by comma or whitespace
            taxon_group = set(re.split(r"[,\s]+", clean_line))
            taxa_groups.append(taxon_group)

    # check 1: definition file must have 1 or 2 taxa groups
    if not taxa_groups:
        raise ValueError(f"Definition file is blank!")
    if len(taxa_groups) > 2:
        raise ValueError(f"Definition file contains more than 2 taxa groups!")
    # if only one group is provided, infer the second group using the input tree
    leaves = set(tree.get_leaf_names())  # assuming leaves are correct as reference
    if len(taxa_groups) == 1:
        taxa_groups.append(leaves - taxa_groups[0])

    # check 2: the two groups of taxa must be non-overlapping, i.e. is there a taxon in both groups?
    if taxa_groups[0] & taxa_groups[1]:
        raise ValueError(f"Defined groups have overlapping items: {taxa_groups[0] & taxa_groups[1]}.")

    # check 3: leaves in the tree must all be in the definition file
    all_taxa = taxa_groups[0] | taxa_groups[1]
    unique_leaves = leaves - all_taxa
    if unique_leaves:
        raise ValueError(f"Some tree leaves do not exist in the definition file: {unique_leaves}.")

    # check 4: all definition file taxa must exist in the tree file
    unique_taxa = all_taxa - leaves
    if unique_taxa:
        raise ValueError(f"Some taxa do not exist in the tree: {unique_taxa}.")

    return taxa_groups

def get_simulation_trees(master_tree, leaf_groups):

    # probe for non-root subtree
    for leaf_group in leaf_groups:
        common_ancestor = master_tree.get_common_ancestor(*leaf_group)

        # returning the potentially rerooted master tree
        # ensures that the first sequences in A_ab alignment correspond
        # to those in A_a alignment, and that the last sequences in
        # A_ab alignment correspond to those in A_b alignment
        # This is necessary for when we merge alignments,
        # the operation of which is sensitive to sequence order
        if not common_ancestor.is_root():
            master_tree.set_outgroup(common_ancestor)
            a_tree, b_tree = master_tree.get_children()
            return master_tree, a_tree, b_tree

if __name__ == "__main__":
    main(args)
