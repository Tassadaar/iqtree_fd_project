#!/usr/bin/env python

import sys
import argparse
import subprocess
import itertools
from ete3 import Tree, TreeNode


def get_ref_subtrees(master_tree, def_address):
    leaf_groups = []
    a_tree = None

    # store definition file in memory
    # as a list of two lists of taxa
    with open(def_address, "r") as def_file:
        # this possible as a list comprehension?
        # suggestion:
        # leaf_groups = [ line.split() for line in def_file ]
        for line in def_file:
            leaf_groups.append(line.split())

    # probe for non-root subtree
    for leaf_group in leaf_groups:
        common_ancestor = master_tree.get_common_ancestor(*leaf_group)

        # how about
        # if not common_ancestor.is_root():
        #    master_tree.set_outgroup(common_ancestor)
        #    a_tree, b_tree = master_tree.get_children()
        #
        # return a_tree, b_tree
        if common_ancestor.is_root():
            continue

        a_tree = common_ancestor

    master_tree.set_outgroup(a_tree)
    b_tree = master_tree.get_children()[1]

    return a_tree, b_tree


# why not simply use BioPython for loading,
# splitting, and writing the alignment files?
def write_alignment_partitions(alignment_address, new_a_tree, new_b_tree):
    alignments = {}

    # parsing fasta file
    with open(alignment_address, "r") as alignment_file:
        current_taxon = ""
        current_alignment = []

        for line in alignment_file:

            if line.startswith(">"):

                if current_taxon and current_alignment:
                    current_alignment[-1] = current_alignment[-1].strip()
                    alignments[current_taxon] = "".join(current_alignment)
                    current_alignment = []

                current_taxon = line[1:-1]
                continue

            current_alignment.append(line)

        if current_taxon and current_alignment:
            alignments[current_taxon] = "".join(current_alignment)

    # write alignment for subtree a
    with open("test_a.aln", "w") as new_a_aln_file:
        a_leaves = new_a_tree.get_leaf_names()

        for taxon, alignment in alignments.items():

            if taxon in a_leaves:
                new_a_aln_file.write(f">{taxon}\n")
                new_a_aln_file.write(f"{alignment}\n")

    # write alignment for subtree b
    with open("test_b.aln", "w") as new_b_aln_file:
        b_leaves = new_b_tree.get_leaf_names()

        for taxon, alignment in alignments.items():

            if taxon in b_leaves:
                new_b_aln_file.write(f">{taxon}\n")
                new_b_aln_file.write(f"{alignment}\n")


def fix_topology(input_tree, reference_tree):
    tree_copy = input_tree.copy("deepcopy")

    # lists of two items, each item a list of taxa
    copied_leaves    = [ child.get_leaf_names() for child in tree_copy.get_children()      ]
    reference_leaves = [ child.get_leaf_names() for child in reference_tree.get_children() ]

    for leaf_group in copied_leaves:

        # check if group of leafs is found as
        # in-group or out-group in reference tree
        ## does list occurr in list of lists?
        if leaf_group not in reference_leaves:
            continue

        # scenario 1: leaf group is in ref tree but is only one taxon
        if len(leaf_group) == 1:
            tree_copy.set_outgroup(leaf_group[0])
            return tree_copy

        # scenario 2: leaf group is in ref tree but is multiple taxa
        # ?? get_common_ancestor() requires multiple arguments ??
        outgroup = tree_copy.get_common_ancestor(*leaf_group)
        tree_copy.set_outgroup(outgroup)
        return tree_copy

    # scenario 3: none of the leaf groups match the reference tree
    # hence we need to recursively climb the subtree until
    # a leaf group is found that matches an in- or out-group
    # of the reference tree
    def find_outgroup(node):
        parent = node.up
        # stop condition: we found a node that matches
        # an in-group or out-group found in reference tree
        if parent.get_leaf_names() in reference_leaves:
            return parent
        else:
            return find_outgroup(parent)

    # make this a two liner to increase readability
    tree_copy.set_outgroup(find_outgroup(tree_copy.get_farthest_leaf()[0]))
    return tree_copy


def run_iqtree_a(tree_file, alignment, prefix, model):

    # meeting note: have consistent iqtree arguments
    # aim for iqtree 2 style arguments??
    iqtree_command = [
        "iqtree",
        "-nt", "2", # meeting note: set number of cpus as command line argument
        "-s", alignment,
        "-te", tree_file,
        "-m", f"LG+{model}+G", # meeting note: make number of G categories user defined on command line, also {model} should include 'LG+' (aka, user defined input)
        "-mwopt",
        "-prec", "10",
        "--prefix", prefix # --prefix does not work with my version of iqtree. Use -pre instead?
    ]

    subprocess.run(iqtree_command)


def get_info(in_file):
    weights = {}
    alpha = 0

    with open(in_file, "r") as file:
        section_found = False

        for line in file:

            if "No  Component      Rate    Weight   Parameters" in line:
                section_found = True
                continue

            if line == "\n":
                section_found = False

            if section_found:
                # a little bit more readable
                category, _, _, weight, _  = line.split()
                weights[category] = float(weight)

            if "Gamma shape alpha:" in line:
                words = line.split()
                alpha = float(words[3])
                break

    return weights, alpha


def get_averages():
    a_weights, a_alpha = get_info("test_subset1.iqtree")
    b_weights, b_alpha = get_info("test_subset2.iqtree")

    # meeting note: means should be the weighted means, weighted by number of taxa
    avg_weights = {category : (a_weight + b_weights[category]) / 2 for category, a_weight in a_weights.items()}
    avg_alpha = (a_alpha + b_alpha) / 2

    return avg_weights, avg_alpha


def write_nexus_file(weights, model):
    out_freqs = []
    out_model = []

    # generate frequency section
    # we take the model frequencies from a iqtree data file
    with open("data/modelmixtureCAT.nex", "r") as models:
        section_found = False

        for line in models:

            if f"CAT-{model} profile mixture model" in line:
                section_found = True
                continue

            if section_found:

                # the section ends if blank line is encountered
                if line == "\n":
                    break

                words = line.split()
                # rename model to fundi_<model>
                words[1] = f"fundi_{words[1]}"
                # store model frequencies in memory
                out_freqs.append(words)

    # generate model section
    # why define these as lists? why not define them as strings?
    # meeting note: poisson+G+ can be removed
    weight_line = [ f"model fundi_{model} = POISSON+G+FMIX{{"    ]
    # meeting note: model_line is not necessary!
    model_line  = [ f"model fundi_{model}Opt = POISSON+G+FMIX{{" ]

    # we take the previously calculated model weights from memory
    ## this code is a little messy
    ## how about
    # 
    # last_category = list( weights.keys() )[-1]
    # for category, weight in weights.items():
    #    if category != last_category:
    #        weight_line.append(f"fundi_{model}pi{category}:1:{weight}" + ",")
    #        model_line.append( f"fundi_{model}pi{cateogry}" + ",")
    #    else:
    #        weight_line.append(f"fundi_{model}pi{category}:1:{weight}" + "};")
    #        model_line.append( f"fundi_{model}pi{cateogry}" + "};")
    #
    # or if you want to have weight_line be a string
    # last_category = list( weights.keys() )[-1]
    # for category, weight in weights.items():
    #    if category != last_category:
    #        weight_line += f"fundi_{model}pi{category}:1:{weight}" + ","
    #        model_line  += f"fundi_{model}pi{cateogry}" + ","
    #    else:
    #        weight_line += f"fundi_{model}pi{category}:1:{weight}" + "};"
    #        model_line  += f"fundi_{model}pi{cateogry}" + "};"
    for index, weight in weights.items():
        weight_line.append(f"fundi_{model}pi{index}:1:{weight}" + ("," if index != list(weights.keys())[-1] else "};"))
        model_line.append(f"fundi_{model}pi{index}" + ("," if index != list(weights.keys())[-1] else "};"))

    out_model.append(weight_line)
    out_model.append(model_line)

    # write to file
    with open("test_nex.nex", "w") as nex_file:
        nex_file.write("#nexus\nbegin models;\n")

        for line in out_freqs:
            nex_file.write(" ".join(line) + "\n")

        for line in out_model:
            # if out_model is a list of strings, instead of lists of lists
            # nex_file.write(line + \n")
            nex_file.write("".join(line) + "\n")

        nex_file.write("end;")


def run_iqtree_b(trees, avg_alpha, model):
    commands = []

    i = 1
    for tree in trees:

        tree.render(f"test_{i}.png")

        # write tree with
        # t.write(format=1, outfile="new_tree.nw")
        with open(f"test_{i}.tree", "w") as tree_file:
            tree_file.write(tree.write())

        # this iqtree command should invoke fundi
        # meeting note: at least iqtree2 v2.2 is necessary
        iqtree_cmd = 
            "iqtree",
            "-s", "data/Dandan/toy.aln", # meeting note: provide alignment as argument on command line
            "--tree-fix", tree_file.name,
            "-m", f"LG+fundi_{model}+G{{{avg_alpha}}}",
            "--mdef", "test_nex.nex",
            "-T", "8", # meeting note: set number of cpus as argument on command line
            "-blfix",
            "--prefix", f"test_{i}",# --prefix does not work with my version of iqtree. Use -pre instead?
            "-prec", "10",
        ]

        commands.append(iqtree_cmd)
        i += 1

    # second iqtree run
    for command in commands:
        subprocess.run(command)


def main(args):

    # check if this kind of try / except block is
    # the recommended way of doing things
    try:
        master_tree = Tree(args.tree)
        master_tree.set_outgroup(master_tree.get_children()[0]) # the master tree needs to be rooted for ete3

        # should it be rooted though?
        assert len(master_tree.get_children()) == 2, f"Master tree must be rooted!\n {master_tree}"

        # just use args.xxxx straight up,
        # rather than creating new memory variables
        full_alignment = args.alignment
        def_file = args.definition
        model = args.mixture_model

        # split master tree into two TreeNode objects
        a_tree, b_tree = get_ref_subtrees(master_tree, def_file)

        # write a newick file in a single line of code with
        # e.g., t.write(format=1, outfile="new_tree.nw")
        with open("test_a.tree", "w") as a_tree_file:
            a_tree_file.write(new_a_tree.write())

        with open("test_b.tree", "w") as b_tree_file:
            b_tree_file.write(new_b_tree.write())

        # generates 'test_a.aln' and 'test_b.aln'
        write_alignment_partitions(full_alignment, new_a_tree, new_b_tree)

        # first iqtree execution
        # perhaps make named arguments for this function?
        # just to increase code clarity
        # e.g. run_iqtree_a(fixed_tree='test_a.tree', alignment='test_a.aln', prefix='test_subset1', model=args.mixture_model)
        run_iqtree_a("test_a.tree", "test_a.aln", "test_subset1", model)
        run_iqtree_a("test_b.tree", "test_b.aln", "test_subset2", model)

        # do we want to use the fix_topology() function here??
        # the iqtree output trees are unrooted, so before we
        # stitch them together we have to make sure they are
        # properly rooted

        # meeting note: take branch lengths from iqtree outputs!!
        a_branch_length = new_a_tree.get_children()[0].dist + new_a_tree.get_children()[1].dist
        b_branch_length = new_b_tree.get_children()[0].dist + new_b_tree.get_children()[1].dist

        print(f"branch a: {a_branch_length}, branch b: {b_branch_length}\n")

        # make a list of increments
        # [0.1, 0.2, 0.3, etc, ...]
        denominator = int(1 / args.increment)
        proportions = [ x / denominator for x in range(1, denominator) ]

        trees = []
        # get alpha and beta from a cartesian product of proportions
        for alpha, beta in itertools.product(proportions, repeat=2):

            # set new branch lengths for a
            new_a_tree.get_children()[0].dist = a_branch * alpha
            new_a_tree.get_children()[1].dist = a_branch * (1 - alpha)

            # set new branch lengths for b
            new_b_tree.get_children()[0].dist = b_branch * beta
            new_b_tree.get_children()[1].dist = b_branch * (1 - beta)

            # stitch subtrees back together into master tree
            new_master_tree = TreeNode(dist=0.1)
            new_master_tree.add_child(new_a_tree.copy("deepcopy"))
            new_master_tree.add_child(new_b_tree.copy("deepcopy"))

            assert new_master_tree.robinson_foulds(master_tree)[0] == 0, "The new master tree is not the same!"

            trees.append(new_master_tree)

        # print number of trees generated
        print(f"{len(trees)} tree were generated")

        # get average weights and alpha
        avg_weights, avg_alpha = get_averages()

        # generate nexus file:
        write_nexus_file(avg_weights, model)

        # second iqtree execution
        run_iqtree_b(trees, avg_alpha, model)

        # meeting note: generate a little summary report, 
        # select the tree with the highest likelihood

    except AssertionError as e:
        print(f"Oops! {e}")

    except NameError as e:
        print(f"Panda-monium! {e}")

    finally:
        print("\nSystem exiting...")
        sys.exit()


# put parser on top, but main function down here
# let main function go
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Tree mapper")

    parser.add_argument("-te", "--tree", required=True, 
                        help="Tree file in newick format, must be rooted")

    parser.add_argument("-s", "--alignment", required=True, 
                        help="Alignment in fasta format")

    # meeting note: add validation function for definition file
    # does the taxa split in def correspond to a bipartition in the master tree?
    # if not kill execution and tell user to fix its definition file
    parser.add_argument("-d", "--definition", required=True,
                        help="Definition file that splits the tree by FunDi branch")

    parser.add_argument("-m", "--mixture_model", required=False, default="C10",
                        help="Mixture model to be used with iqtree")

    parser.add_argument("-i", "--increment", required=False, default=0.1,
                        help="Metric to control branch length variance, default is 0.1")

    # # emulating commandline arguments for development
    # sys.argv = [
    #     "Mapper.py",
    #     "-te", "data/Dandan/rooted_toy.newick",
    #     "-d", "data/Dandan/toy.def",
    #     "-s", "data/Dandan/toy.aln",
    # ]

    arguments = parser.parse_args()
    main(arguments)
