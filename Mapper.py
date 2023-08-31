import sys
import argparse
import subprocess
from itertools import product
from ete3 import Tree, TreeNode


def get_ref_subtrees(master_tree, def_address):
    leaf_groups = []
    a_tree = None

    with open(def_address, "r") as def_file:

        for line in def_file:
            leaf_groups.append(line.split())

    # probe for non-root subtree
    for leaf_group in leaf_groups:
        tree = master_tree.get_common_ancestor(*leaf_group)

        if tree.is_root():
            continue

        a_tree = tree

    master_tree.set_outgroup(a_tree)
    b_tree = master_tree.get_children()[1]

    return a_tree, b_tree


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
    copied_leaves = [child.get_leaf_names() for child in tree_copy.get_children()]
    reference_leaves = [child.get_leaf_names() for child in reference_tree.get_children()]

    for leaf_group in copied_leaves:

        if leaf_group not in reference_leaves:
            continue

        # scenario 1: leaf group is valid but only one taxon
        if len(leaf_group) == 1:
            tree_copy.set_outgroup(leaf_group[0])
            return tree_copy

        # scenario 2: leaf group is valid but multiple taxa
        outgroup = tree_copy.get_common_ancestor(*leaf_group)
        tree_copy.set_outgroup(outgroup)
        return tree_copy

    # scenario 3: none of the leaf groups of root is valid
    def find_outgroup(node):
        parent = node.up

        if parent.get_leaf_names() in reference_leaves:
            return parent
        else:
            return find_outgroup(parent)

    tree_copy.set_outgroup(find_outgroup(tree_copy.get_farthest_leaf()[0]))
    return tree_copy


def run_iqtree_a(tree_file, alignment, prefix, model):

    iqtree_command = [
        "iqtree",
        "-nt", "2",
        "-s", alignment,
        "-te", tree_file,
        "-m", f"LG+{model}+G",
        "-mwopt",
        "-prec", "10",
        "--prefix", prefix
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
                words = line.split()
                weights[words[0]] = float(words[3])

            if "Gamma shape alpha:" in line:
                words = line.split()
                alpha = float(words[3])
                break

    return weights, alpha


def get_averages():
    a_log = "test_subset1.iqtree"
    b_log = "test_subset2.iqtree"

    a_weights, a_alpha = get_info(a_log)
    b_weights, b_alpha = get_info(b_log)

    avg_weights = {key: (a_weight + b_weights[key]) / 2 for key, a_weight in a_weights.items()}
    avg_alpha = (a_alpha + b_alpha) / 2

    return avg_weights, avg_alpha


def write_nexus_file(weights, model):
    out_freqs = []
    out_model = []

    # generate frequency section
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
                words[1] = f"fundi_{words[1]}"
                out_freqs.append(words)

    # generate model section
    weight_line = [f"model fundi_{model} = POISSON+G+FMIX{{"]
    model_line = [f"model fundi_{model}Opt = POISSON+G+FMIX{{"]

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
            nex_file.write("".join(line) + "\n")

        nex_file.write("end;")


def run_iqtree_b(trees, avg_alpha, model):
    commands = []

    i = 1
    for tree in trees:

        tree.render(f"test_{i}.png")

        with open(f"test_{i}.tree", "w") as tree_file:
            tree_file.write(tree.write())

        iqtree_cmd = [
            "iqtree",
            "-s", "data/Dandan/toy.aln",
            "--tree-fix", tree_file.name,
            "-m", f"LG+fundi_{model}+G{{{avg_alpha}}}",
            "--mdef", "test_nex.nex",
            "-T", "8",
            "-blfix",
            "--prefix", f"test_{i}",
            "-prec", "10",
        ]

        commands.append(iqtree_cmd)
        i += 1

    # second iqtree run
    for command in commands:
        subprocess.run(command)


def main(args):

    try:
        master_tree = Tree(args.tree)
        master_tree.set_outgroup(master_tree.get_children()[0]) # the master tree needs to be rooted for ete3

        assert len(master_tree.get_children()) == 2, f"Master tree must be rooted!\n {master_tree}"

        full_alignment = args.alignment
        def_file = args.definition
        model = args.mixture_model
        denominator = int(1 / args.increment)
        a_tree, b_tree = get_ref_subtrees(master_tree, def_file)

        # ------------------------------------section below is largely deprecated---------------------------------------
        # check if a_tree and b_tree have the same topology
        # subtrees used to be user input during development
        # "Do you think this is worth keeping? I doubt it."
        for subtree in master_tree.get_children():

            if subtree.get_leaf_names() == a_tree.get_leaf_names():
                sub_a_tree = subtree

            if subtree.get_leaf_names() == b_tree.get_leaf_names():
                sub_b_tree = subtree

        if a_tree.robinson_foulds(sub_a_tree)[0] != 0:
            new_a_tree = fix_topology(a_tree, sub_a_tree)
        else:
            new_a_tree = a_tree

        if b_tree.robinson_foulds(sub_b_tree)[0] != 0:
            new_b_tree = fix_topology(b_tree, sub_b_tree)
        else:
            new_b_tree = b_tree
        # --------------------------------------------------------------------------------------------------------------

        with open("test_a.tree", "w") as a_tree_file:
            a_tree_file.write(new_a_tree.write())

        with open("test_b.tree", "w") as b_tree_file:
            b_tree_file.write(new_b_tree.write())

        write_alignment_partitions(full_alignment, new_a_tree, new_b_tree)

        # first iqtree execution
        run_iqtree_a("test_a.tree", "test_a.aln", "test_subset1", model)
        run_iqtree_a("test_b.tree", "test_b.aln", "test_subset2", model)

        trees = []
        proportions = [x / denominator for x in range(1, denominator)]
        a_branch = new_a_tree.get_children()[0].dist + new_a_tree.get_children()[1].dist
        b_branch = new_b_tree.get_children()[0].dist + new_b_tree.get_children()[1].dist

        print(f"branch a: {a_branch}, branch b: {b_branch}\n")

        # get alpha and beta from a cartesian product of proportions
        for alpha, beta in product(proportions, repeat=2):

            # set new branch lengths for a
            new_a_tree.get_children()[0].dist = a_branch * alpha
            new_a_tree.get_children()[1].dist = a_branch * (1 - alpha)

            # set new branch lengths for b
            new_b_tree.get_children()[0].dist = b_branch * beta
            new_b_tree.get_children()[1].dist = b_branch * (1 - beta)

            # reconstruct master tree
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

    except AssertionError as e:
        print(f"Oops! {e}")

    except NameError as e:
        print(f"Panda-monium! {e}")

    finally:
        print("\nSystem exiting...")
        sys.exit()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Tree mapper")

    parser.add_argument("-te", "--tree", required=True, help="Tree file in newick format, must be rooted")
    parser.add_argument("-s", "--alignment", required=True, help="Alignment in fasta format")
    parser.add_argument("-d", "--definition", required=True,
                        help="Definition file that splits the tree by FunDi branch")
    parser.add_argument("-m", "--mixture_model", required=False, default="C10",
                        help="Mixture model to be used with iqtree")
    parser.add_argument("-i", "--increment", required=False, default=0.1,
                        help="Metric to control branch length variance, default is 0.1")

    # emulating commandline arguments for development
    sys.argv = [
        "Mapper.py",
        "-te", "data/Dandan/gstoy.newick",
        "-d", "data/Dandan/toy.def",
        "-s", "data/Dandan/toy.aln",
    ]

    arguments = parser.parse_args()
    main(arguments)
