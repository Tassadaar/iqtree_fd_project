import os
import sys
import argparse
import subprocess
from itertools import product
from ete3 import Tree, TreeNode
from Bio import AlignIO


def validate_alignment(tree, alignment_address):
    # parsing fasta file
    alignment = AlignIO.read(alignment_address, "fasta")
    seq_ids = set(record.id for record in alignment)

    for taxon in tree.get_leaf_names():

        if taxon not in seq_ids:
            raise ValueError(f"Taxon {taxon} does not have alignment information.")

    return alignment


def validate_def_file(tree, def_address):
    # store definition file in memory
    # as a list of two lists of taxa
    with open(def_address, "r") as def_file:
        leaf_groups = [line.split() for line in def_file]

    if len(leaf_groups) != 2:
        raise ValueError(f"Definition file is not in the right format.")

    leaves = leaf_groups[0] + leaf_groups[1]
    remainder_leaves = leaves.copy()

    for taxon in tree.get_leaf_names():

        if taxon not in leaves:
            raise ValueError(f"Taxon {taxon} is undefined in the definition file")

        remainder_leaves.remove(taxon)

    if remainder_leaves:

        if len(remainder_leaves) == 1:
            raise ValueError(f"Taxon {{{remainder_leaves.pop()}}} in the definition file do not exist in the tree.")
        else:
            raise ValueError(f"Taxa {{{', '.join(remainder_leaves)}}} in the definition file do not exist in the tree.")

    return leaf_groups


def get_ref_subtrees(master_tree, leaf_groups):
    # this is good practice, but it currently breaks the code
    # master_tree_copy = master_tree.copy("deepcopy")
    a_tree = None
    b_tree = None

    # probe for non-root subtree
    for leaf_group in leaf_groups:
        common_ancestor = master_tree.get_common_ancestor(*leaf_group)

        if not common_ancestor.is_root():
            master_tree.set_outgroup(common_ancestor)
            a_tree, b_tree = master_tree.get_children()

    return a_tree, b_tree


def write_alignment_partitions(alignment, a_tree, b_tree):
    a_leaves = a_tree.get_leaf_names()
    b_leaves = b_tree.get_leaf_names()

    # write alignment for subtree a
    with open("test_a.aln", "w") as a_aln_file:
        
        for record in alignment:

            if record.id in a_leaves:
                a_aln_file.write(f">{record.id}\n")
                a_aln_file.write(f"{record.seq}\n")

    # write alignment for subtree b
    with open("test_b.aln", "w") as b_aln_file:

        for record in alignment:

            if record.id in b_leaves:
                b_aln_file.write(f">{record.id}\n")
                b_aln_file.write(f"{record.seq}\n")


def run_iqtree_a(tree_file, alignment_address, prefix, model):

    iqtree_command = [
        "iqtree",
        "-nt", "2",
        "-s", alignment_address,
        "-te", tree_file,
        "-m", f"LG+{model}+G",
        "-mwopt",
        "-prec", "10",
        "--prefix", prefix
    ]

    subprocess.run(iqtree_command)


def validate_iqtree_generated_tree(iqtree_file, ref_tree):

    def fix_topology(input_tree):
        ref_outgroup_leaves = ref_tree.get_children()[0].get_leaf_names()

        if len(ref_outgroup_leaves) == 1:
            outgroup = ref_outgroup_leaves[0]
        else:
            outgroup = input_tree.get_common_ancestor(*ref_outgroup_leaves)

            if outgroup.is_root():

                for leaf in input_tree.get_leaf_names():

                    if leaf not in ref_outgroup_leaves:
                        input_tree.set_outgroup(leaf)
                        outgroup = input_tree.get_common_ancestor(*ref_outgroup_leaves)

        input_tree.set_outgroup(outgroup)

    iqtree = None

    with open(iqtree_file, "r") as file:
        section_found = False

        for line in file:

            if "Tree in newick format:" in line:
                section_found = True
                continue

            if line == "\n":
                continue

            if section_found:
                iqtree = Tree(line)

                break

    # root iqtree randomly
    iqtree.set_outgroup(iqtree.get_children()[0])

    if iqtree.robinson_foulds(ref_tree)[0] != 0:
        fix_topology(iqtree)

    return iqtree


def calculate_weights(a_tree, b_tree):
    a_taxa_count = len(a_tree.get_leaf_names())
    b_taxa_count = len(b_tree.get_leaf_names())
    total_taxa_count = a_taxa_count + b_taxa_count

    return a_taxa_count / total_taxa_count, b_taxa_count / total_taxa_count


def calculate_weighted_average_mixture_weights(a_log_file, b_log_file, a_weight, b_weight):

    # get weights from iqtree log file, returns empty set if not found
    def get_mixture_weights(iqtree_file):
        mixture_weights = {}

        if not os.path.isfile(iqtree_file):
            raise FileNotFoundError(f"'{iqtree_file}' does not exist!")

        with open(iqtree_file, "r") as file:
            section_found = False

            for line in file:

                if "No  Component      Rate    Weight   Parameters" in line:
                    section_found = True
                    continue

                if line == "\n":
                    section_found = False

                if section_found:
                    words = line.split()
                    mixture_weights[words[0]] = float(words[3])

                if "Gamma shape alpha:" in line:
                    break

        return mixture_weights

    a_mixture_weights = get_mixture_weights(a_log_file)
    b_mixture_weights = get_mixture_weights(b_log_file)

    if not a_mixture_weights:
        raise ValueError(f"Cannot extract weights, check if '{a_log_file}' is formatted correctly!")

    if not b_mixture_weights:
        raise ValueError(f"Cannot extract weights, check if '{a_log_file}' is formatted correctly!")

    avg_mixture_weights = {
        key: (a_mixture_weight * a_weight + b_mixture_weights[key] * b_weight) / 2
        for key, a_mixture_weight in a_mixture_weights.items()
    }

    return avg_mixture_weights


def calculate_weighted_average_alpha(a_log_file, b_log_file, a_weight, b_weight):

    # get alpha from iqtree log file, returns None if not found
    def get_alpha(iqtree_file):
        alpha = None

        if not os.path.isfile(iqtree_file):
            raise FileNotFoundError(f"'{iqtree_file}' does not exist!")

        with open(iqtree_file, "r") as file:

            for line in file:

                if "Gamma shape alpha:" in line:
                    words = line.split()
                    alpha = float(words[3])

                    return alpha

        return alpha

    a_alpha = get_alpha(a_log_file)
    b_alpha = get_alpha(b_log_file)

    if not a_alpha:
        raise ValueError(f"Cannot extract weights, check if '{a_log_file}' is formatted correctly!")

    if not b_alpha:
        raise ValueError(f"Cannot extract weights, check if '{a_log_file}' is formatted correctly!")

    return (a_alpha * a_weight + b_alpha * b_weight) / 2


def write_nexus_file(weights, model):
    file_name = "test_nex.nex"
    out_freqs = []

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

    for index, weight in weights.items():
        weight_line.append(f"fundi_{model}pi{index}:1:{weight}" + ("," if index != list(weights.keys())[-1] else "};"))

    # write to file
    with open(file_name, "w") as nex_file:
        nex_file.write("#nexus\nbegin models;\n")

        for line in out_freqs:
            nex_file.write(" ".join(line) + "\n")

        nex_file.write("".join(weight_line) + "\n")

        nex_file.write("end;")

    return file_name


def run_iqtree_b(trees, alignment_address, avg_alpha, model, nexus_file):
    commands = []

    i = 1
    for tree in trees:

        tree.render(f"test_{i}.png")

        with open(f"test_{i}.tree", "w") as tree_file:
            tree_file.write(tree.write())

        iqtree_cmd = [
            "iqtree",
            "-s", alignment_address,
            "--tree-fix", tree_file.name,
            "-m", f"LG+fundi_{model}+G{{{avg_alpha}}}",
            "--mdef", nexus_file,
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


def generate_summary(tree_count):
    likelihoods = {}

    i = 1
    while i <= tree_count:

        with open(f"test_{i}.iqtree", "r") as iqtree_file:

            for line in iqtree_file:

                if "Log-likelihood of the tree:" not in line:
                    continue

                words = line.split()
                likelihoods[i] = float(words[4])
                continue

        i += 1

    with open("test_summary.txt", "w") as summary_file:

        best_tree = max(likelihoods, key=likelihoods.get)
        summary_file.write(f"Tree {best_tree} has the largest log-likelihood of {likelihoods[best_tree]}.\n\n")

        with open(f"test_{best_tree}.iqtree", "r") as tree_file:
            section_found = False

            for line in tree_file:

                if "NOTE: Tree is UNROOTED" in line:
                    section_found = True

                if "Tree in newick format:" in line:
                    break

                if section_found is True:
                    summary_file.write(line)

        for tree, likelihood in likelihoods.items():
            summary_file.write(f"Log-likelihood of the tree {tree}: {likelihood}\n")


def main(args):

    try:
        master_tree = Tree(args.tree)
        master_tree.set_outgroup(master_tree.get_children()[0])  # the master tree needs to be rooted for ete3

        assert len(master_tree.get_children()) == 2, f"Master tree must be rooted!\n {master_tree}"

        alignment_address = args.alignment
        alignment = validate_alignment(master_tree, alignment_address)
        defined_groups = validate_def_file(master_tree, args.definition)
        model = args.mixture_model
        nexus_address = args.nexus
        denominator = int(1 / args.increment)
        a_tree, b_tree = get_ref_subtrees(master_tree, defined_groups)

        with open("test_a.tree", "w") as a_tree_file:
            a_tree_file.write(a_tree.write())

        with open("test_b.tree", "w") as b_tree_file:
            b_tree_file.write(b_tree.write())

        write_alignment_partitions(alignment, a_tree, b_tree)

        # first iqtree execution
        run_iqtree_a("test_a.tree", "test_a.aln", "test_subset1", model)
        run_iqtree_a("test_b.tree", "test_b.aln", "test_subset2", model)
        a_iqtree_file = "test_subset1.iqtree"
        b_iqtree_file = "test_subset2.iqtree"

        # check if the new trees generated by iqtree have the same topology
        a_tree = validate_iqtree_generated_tree(a_iqtree_file, a_tree)
        b_tree = validate_iqtree_generated_tree(b_iqtree_file, b_tree)

        trees = []
        proportions = [x / denominator for x in range(1, denominator)]
        a_branch = a_tree.get_children()[0].dist + a_tree.get_children()[1].dist
        b_branch = b_tree.get_children()[0].dist + b_tree.get_children()[1].dist

        print(f"branch a: {a_branch}, branch b: {b_branch}\n")

        # get alpha and beta from a cartesian product of proportions
        for alpha, beta in product(proportions, repeat=2):

            # set new branch lengths for a
            a_tree.get_children()[0].dist = a_branch * alpha
            a_tree.get_children()[1].dist = a_branch * (1 - alpha)

            # set new branch lengths for b
            b_tree.get_children()[0].dist = b_branch * beta
            b_tree.get_children()[1].dist = b_branch * (1 - beta)

            # reconstruct master tree
            new_master_tree = TreeNode(dist=0.1)
            new_master_tree.add_child(a_tree.copy("deepcopy"))
            new_master_tree.add_child(b_tree.copy("deepcopy"))

            assert new_master_tree.robinson_foulds(master_tree)[0] == 0, "The new master tree is not the same!"

            trees.append(new_master_tree)

        # print number of trees generated
        print(f"{len(trees)} tree were generated")

        a_weight, b_weight = calculate_weights(a_tree, b_tree)
        avg_alpha = calculate_weighted_average_alpha(a_iqtree_file, b_iqtree_file, a_weight, b_weight)

        if not nexus_address:
            avg_mixture_weights = calculate_weighted_average_mixture_weights(a_iqtree_file, b_iqtree_file, a_weight, b_weight)
            # generate nexus file:
            nexus_address = write_nexus_file(avg_mixture_weights, model)

        # second iqtree execution
        run_iqtree_b(trees, alignment_address, avg_alpha, model, nexus_address)

        # To get the number of trees generated, we take number of proportions to the power of 2
        generate_summary(len(proportions) ** 2)

    except AssertionError as e:
        print(f"Oops! {e}")

    except NameError as e:
        print(f"Panda-monium! {e}")

    except ValueError as e:
        print(f"Fatal: {e} Check if inputs are valid!")

    except FileNotFoundError as e:
        print(f"File not found! {e}")

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
    parser.add_argument("-mdef", "--nexus", required=False, default=None,
                        help="Nexus file to be used with iqtree")

    # emulating commandline arguments for development
    sys.argv = [
        "Mapper.py",
        "-te", "data/Dandan/toy.newick",
        "-d", "data/Dandan/toy.def",
        "-s", "data/Dandan/toy.aln",
    ]

    arguments = parser.parse_args()
    main(arguments)
    # tree = Tree("data/Dandan/toy.subset1.newick")
    # new_tree = validate_iqtree_generated_tree("data/Dandan/toy.subset1.aln.iqtree", tree)
    # print(new_tree)