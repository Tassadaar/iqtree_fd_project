import os
import sys
import argparse
import subprocess
from itertools import product
from ete3 import Tree, TreeStyle, TreeNode
from Bio import AlignIO


def validate_alignment(tree, alignment_file):
    """
    Does the alignment have the exact same set of taxa
    as the tree?
    """
    # parsing fasta file
    alignment = AlignIO.read(alignment_file, "fasta")
    seq_ids = set(record.id for record in alignment)
    seen_taxa = set()

    for taxon in tree.get_leaf_names():

        if taxon not in seq_ids:
            raise ValueError(f"Tree taxon {taxon} does not exist in the alignment file")

        if taxon in seen_taxa:
            raise ValueError(f"Tree taxon {taxon} is duplicated in the alignment file")

        seen_taxa.add(taxon)

    return alignment


def validate_def_file(tree, def_file):
    # store definition file in memory
    # as a list of two sets of taxa
    # this ensures uniqueness, assuming that order does not matter with def files
    with open(def_file, "r") as file:
        taxa_groups = [set(line.split()) for line in file]

    # check 1: definition file must have exactly two groups of taxa
    if len(taxa_groups) != 2:
        raise ValueError(f"Definition file does not have exactly two groups of taxa")

    # check 2: the two groups of taxa must be non-overlapping, i.e. is there a taxon in both groups?
    set_intercept = taxa_groups[0] & taxa_groups[1]

    if set_intercept:
        raise ValueError(f"Defined groups have overlapping items: {set_intercept}.")

    # check 3: leaves in the tree must all be in the definition file
    all_taxa = taxa_groups[0] | taxa_groups[1]
    leaves = set(tree.get_leaf_names())
    unique_leaves = leaves - all_taxa

    if unique_leaves:
        raise ValueError(f"Some tree leaves do not exist in the definition file: {unique_leaves}.")

    # check 4: are there any taxa in the definition file that do not exist in the tree file?
    unique_taxa = all_taxa - leaves

    if unique_taxa:
        raise ValueError(f"Some taxa do not exist in the tree: {unique_taxa}.")

    return taxa_groups


def get_ref_subtrees(master_tree, leaf_groups):
    # this is good practice, but it currently breaks the code
    # master_tree_copy = master_tree.copy("deepcopy")

    # we declare these as None, so that if
    # in the code below they are not replaced with actual objects,
    # we can catch the error later
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
    # TODO: use a_tree.get_leaf_names() directly, rather than storing it seperately
    #  NOTE: we don't want to call the method every loop, so no
    a_leaves = a_tree.get_leaf_names()
    b_leaves = b_tree.get_leaf_names()
    a_alignment = AlignIO.MultipleSeqAlignment([])
    b_alignment = AlignIO.MultipleSeqAlignment([])

    for record in alignment:

        if record.id in a_leaves:
            a_alignment.append(record)

        if record.id in b_leaves:
            b_alignment.append(record)

    AlignIO.write(a_alignment, "test_a.aln", "fasta")
    AlignIO.write(b_alignment, "test_b.aln", "fasta")


def conform_iqtree_tree(iqtree_file, ref_tree):

    # extracting subtree from iqtree file
    iqtree_tree = None
    with open(iqtree_file, "r") as file:
        # section_found = False

        for line in file:

            if "Tree in newick format:" in line:
                next(file)
                iqtree_tree = Tree(next(file))
                # section_found = True
                # continue

            # if line == "\n":
                # continue

            # if section_found:
                # iqtree_tree = Tree(line)
                # break

    # rooting subtree using reference subtree as template
    ref_outgroup_leaves = ref_tree.get_children()[0].get_leaf_names()

    if len(ref_outgroup_leaves) == 1:
        outgroup = ref_outgroup_leaves[0]
    else:
        outgroup = iqtree_tree.get_common_ancestor(*ref_outgroup_leaves)

        if outgroup.is_root():

            for leaf in iqtree_tree.get_leaf_names():

                if leaf not in ref_outgroup_leaves:
                    iqtree_tree.set_outgroup(leaf)
                    outgroup = iqtree_tree.get_common_ancestor(*ref_outgroup_leaves)

    iqtree_tree.set_outgroup(outgroup)


    # # root iqtree randomly
    # iqtree_tree.set_outgroup(iqtree_tree.get_children()[0])

    # if iqtree.robinson_foulds(ref_tree)[0] != 0:
    # fix_topology(iqtree_tree)

    return iqtree_tree


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
        key: a_mixture_weight * a_weight + b_mixture_weights[key] * b_weight
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

    return a_alpha * a_weight + b_alpha * b_weight


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
    weight_line = f"model fundi_{model} = FMIX{{"
    last_category = list(weights.keys())[-1]

    for category, weight in weights.items():

        if category != last_category:
            weight_line += f"fundi_{model}pi{category}:1:{weight},"
        else:
            weight_line += f"fundi_{model}pi{category}:1:{weight}" + "};"

    # write to file
    with open(file_name, "w") as nex_file:
        nex_file.write("#nexus\nbegin models;\n")

        for line in out_freqs:
            nex_file.write(" ".join(line) + "\n")

        nex_file.write(weight_line + "\n")

        nex_file.write("end;")

    return file_name


def run_iqtree_b(trees, alignment_address, avg_alpha, model, nexus_file, cores, leaves):
    commands = []

    i = 1
    for tree in trees:

        tree.render(f"test_{i}.png")

        with open(f"test_{i}.tree", "w") as tree_file:
            tree_file.write(tree.write())

        iqtree_cmd = [
            "iqtree2",
            "-s", alignment_address,
            "--tree-fix", tree_file.name,
            "-m", f"LG+fundi_{model}+G{{{avg_alpha}}}",
            "--mdef", nexus_file,
            "-nt", cores,
            "--prefix", f"test_{i}",
            "-prec", "10",
            "-blfix",
            "--fundi", f"{','.join(leaves)},estimate",
            "-redo",
            "--quiet"
        ]

        commands.append(iqtree_cmd)
        i += 1

    # second iqtree run
    for command in commands:
        print(f"Running iqtree funDi on Tree {commands.index(command) + 1} out of {len(commands)}...")
        subprocess.run(command)


def generate_summary(tree_count):

    # get attributes
    attributes = {}

    i = 1
    while i <= tree_count:
        attribute = [0, 0, 0]  # a list of two attributes, fundi log-likelihood, rho value and central branch length

        with open(f"test_{i}.log", "r") as iqtree_file:

            for line in iqtree_file:

                if "Best FunDi parameter rho:" in line:
                    words = line.split()
                    attribute[1] = float(words[4])
                    continue

                if "Best FunDi central branch length:" in line:
                    words = line.split()
                    attribute[2] = float(words[5])
                    continue

                if "FunDi log-likelihood:" in line:
                    words = line.split()
                    attribute[0] = float(words[2])
                    break

        attributes[i] = attribute
        i += 1

    # get the best tree based on funDi log-likelihood
    best_tree_index = max(attributes, key=lambda key: attributes[key][0])
    best_tree = None

    with open(f"test_{best_tree_index}.iqtree", "r") as tree_file:
        section_found = False

        for line in tree_file:

            if "Tree in newick format:" in line:
                section_found = True
                continue

            if section_found is True:

                if line != "\n":
                    best_tree = Tree(line)
                    break

    # style the tree
    tree_style = TreeStyle()
    tree_style.show_leaf_name = True
    tree_style.show_branch_length = True
    tree_style.branch_vertical_margin = 10

    # give each internal node an explicit name
    for node in best_tree.traverse():

        if not node.is_leaf():
            node.name = "node"

    best_tree.render(file_name=f"test_{best_tree_index}.png", tree_style=tree_style, units="px", w=800, h=1000)

    # print to summary file
    print("Generating summary...")

    with open("test_summary.txt", "w") as summary_file:

        summary_file.write(
            f"Tree {best_tree_index} has the largest funDi log-likelihood of {attributes[best_tree_index][0]}.\n"
            f"rho: {attributes[best_tree_index][1]}.\n"
            f"Central branch length: {attributes[best_tree_index][2]}\n"
        )

        # print tree with branch lengths
        summary_file.write(f"{best_tree.get_ascii(attributes=['name', 'dist'], show_internal=True)}\n\n")
        summary_file.write(f"See \"test_{best_tree_index}.png\" for a tree illustration.\n\n")

        for tree, attribute in attributes.items():
            summary_file.write(f"funDi Log-likelihood of the tree {tree}: {attribute[0]}; "
                               f"rho: {attribute[1]}; "
                               f"central branch length: {attribute[2]}\n")

    print("Summary generated under 'test_summary.txt'.")


def main(args):

    try:
        master_tree = Tree(args.tree)

        # NOTE: is this absolutely necessary?
        # master_tree.set_outgroup(master_tree.get_children()[0])  # the master tree needs to be rooted for ete3
        # assert len(master_tree.get_children()) == 2, f"Master tree must be rooted!\n {master_tree}"

        # TODO: remove "redundant" variable names for arg arguments
        alignment_address = args.alignment
        alignment = validate_alignment(master_tree, alignment_address)
        defined_groups = validate_def_file(master_tree, args.definition)
        model = args.mixture_model
        nexus_address = args.nexus
        cores = args.cores
        a_tree, b_tree = get_ref_subtrees(master_tree, defined_groups)
        a_leaves = a_tree.get_leaf_names()

        # write subtrees into newick files
        ## TODO: use
        ## t.write(format=1, outfile="new_tree.nw")
        ## instead
        with open("test_a.tree", "w") as a_tree_file:
            a_tree_file.write(a_tree.write())

        with open("test_b.tree", "w") as b_tree_file:
            b_tree_file.write(b_tree.write())

        # split alignment into two subalignments
        write_alignment_partitions(alignment, a_tree, b_tree)

        # first iqtree execution
        for subset in ['test_a','test_b']:

            print(f"Running iqtree on subtree {subset}...")
            iqtree_command = [
                "iqtree2",
                "-nt", cores,
                "-s", f"{subset}.aln",
                "-te", f"{subset}.tree",
                "-m", f"LG+{model}+G",
                "-mwopt",
                "-prec", "10",
                "--prefix", f"{subset}",
                "--quiet"
            ]
            subprocess.run(iqtree_command)

           
        # # first iqtree execution
        # print("Running iqtree on subtree a...")
        # run_iqtree_a("test_a.tree", "test_a.aln", "test_subset1", model, cores)
        # print("Running iqtree on subtree b...")
        # run_iqtree_a("test_b.tree", "test_b.aln", "test_subset2", model, cores)
        # a_iqtree_file = "test_subset1.iqtree"
        # b_iqtree_file = "test_subset2.iqtree"

        # check if the new trees generated by iqtree have the same topology
        a_tree = conform_iqtree_tree("test_a.iqtree", a_tree)
        b_tree = conform_iqtree_tree("test_b.iqtree", b_tree)

        trees = []
        denominator = int(1 / float(args.increment))
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
        print(f"{len(trees)} tree were generated\n")

        a_weight, b_weight = calculate_weights(a_tree, b_tree)
        avg_alpha = calculate_weighted_average_alpha("test_a.iqtree", "test_b.iqtree", a_weight, b_weight)

        if not nexus_address:
            avg_mixture_weights = calculate_weighted_average_mixture_weights(
                "test_a.iqtree", "test_b.iqtree", a_weight, b_weight
            )
            # generate nexus file:
            nexus_address = write_nexus_file(avg_mixture_weights, model)

        # second iqtree execution
        run_iqtree_b(trees, alignment_address, avg_alpha, model, nexus_address, cores, a_leaves)

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

    parser.add_argument("-i", "--increment", required=False, default="0.1",
                        help="Metric to control branch length variance, default is 0.1")

    parser.add_argument("-mdef", "--nexus", required=False, default=None,
                        help="Nexus file to be used with iqtree")

    parser.add_argument("-nt", "--cores", required=False, default="2",
                        help="Number of CPU cores to use")

    # emulating commandline arguments for development
    sys.argv = [
        "Mapper.py",
        "-te", "data/Hector/TAB",
        "-d", "data/Hector/def",
        "-s", "data/Hector/conAB1rho60.fa",
        "-i", "0.1"
    ]

    arguments = parser.parse_args()

    main(arguments)
