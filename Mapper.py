import os
import sys
import argparse
import subprocess
import itertools
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
    set_intersect = taxa_groups[0] & taxa_groups[1]

    if set_intersect:
        raise ValueError(f"Defined groups have overlapping items: {set_intersect}.")

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
    # extracting subtree 
    iqtree_tree = Tree(iqtree_file)

    # rooting subtree using reference subtree as template
    ref_outgroup_leaves = ref_tree.get_children()[0].get_leaf_names()

    # if outgroup is a single taxon
    if len(ref_outgroup_leaves) == 1:
        iqtree_tree.set_outgroup(ref_outgroup_leaves[0])
        return iqtree_tree

    # if outgroup is multiple taxa
    outgroup_lca = iqtree_tree.get_common_ancestor(*ref_outgroup_leaves)

    # if the lca of the ref tree outgroup taxa is
    # the root of the new iqtree tree, 
    # the ref outgroup taxa are not monophyletic
    if outgroup_lca.is_root():

        for leaf in iqtree_tree.get_leaf_names():

            # we do an initial reroot here,
            # to ensure that the ref outgroup taxa are monophyletic
            if leaf not in ref_outgroup_leaves:
                iqtree_tree.set_outgroup(leaf)
                break

    # now that they are monophyletic,
    # we can reroot with the reference outgroup taxa
    new_outgroup_lca = iqtree_tree.get_common_ancestor(*ref_outgroup_leaves)
    iqtree_tree.set_outgroup(new_outgroup_lca)

    return iqtree_tree


def calculate_weights(a_tree, b_tree):
    a_taxa_count = len(a_tree.get_leaf_names())
    b_taxa_count = len(b_tree.get_leaf_names())
    total_taxa_count = a_taxa_count + b_taxa_count

    return a_taxa_count / total_taxa_count, b_taxa_count / total_taxa_count


def calculate_weighted_average_mixture_weights(a_iqtree_file, b_iqtree_file, a_weight, b_weight):

    # get weights from iqtree log file, returns empty set if not found
    def get_mixture_weights(iqtree_file):
        mixture_weights = {}

        if not os.path.isfile(iqtree_file):
            raise FileNotFoundError(f"'{iqtree_file}' does not exist!")

        with open(iqtree_file, "r") as file:

            for line in file:

                if "No  Component      Rate    Weight   Parameters" in line:
                    next_line = next(file)

                    while next_line != "\n":
                        words = next_line.split()
                        mixture_weights[words[0]] = float(words[3])
                        next_line = next(file)

                    break

        return mixture_weights

    a_mixture_weights = get_mixture_weights(a_iqtree_file)
    b_mixture_weights = get_mixture_weights(b_iqtree_file)

    # usefulness of these safety nets questionable,
    # but doesn't hurt
    if not a_mixture_weights:
        raise ValueError(f"Cannot extract weights, check if '{a_iqtree_file}' is formatted correctly!")

    if not b_mixture_weights:
        raise ValueError(f"Cannot extract weights, check if '{a_iqtree_file}' is formatted correctly!")

    # mixture_weight = weight of the mixture model class
    # weight = weight used to calculate the weighted average
    avg_mixture_weights = {}

    # a_mixture_weights and b_mixture_weights have the same keys
    for key, a_mixture_weight in a_mixture_weights.items():
        weighted_avg_weight = a_mixture_weight * a_weight + b_mixture_weights[key] * b_weight
        avg_mixture_weights[key] = weighted_avg_weight

    return avg_mixture_weights


def calculate_weighted_average_alpha(a_iqtree_file, b_iqtree_file, a_weight, b_weight):

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

    a_alpha = get_alpha(a_iqtree_file)
    b_alpha = get_alpha(b_iqtree_file)

    if not a_alpha:
        raise ValueError(f"Cannot extract weights, check if '{a_iqtree_file}' is formatted correctly!")

    if not b_alpha:
        raise ValueError(f"Cannot extract weights, check if '{a_iqtree_file}' is formatted correctly!")

    return a_alpha * a_weight + b_alpha * b_weight


def write_nexus_file(weights, model):
    out_freqs = []

    # generate frequency section
    with open("data/modelmixtureCAT.nex", "r") as models:

        for line in models:

            if f"CAT-{model} profile mixture model" in line:
                next_line = next(models)

                while next_line != "\n":
                    words = next_line.split()
                    words[1] = f"fundi_{words[1]}"
                    out_freqs.append(words)
                    next_line = next(models)

                break

    # generate model section
    weight_statements = []

    for category, weight in weights.items():
        weight_statements.append(f"fundi_{model}pi{category}:1:{weight}")

    weight_line = f"model fundi_{model} = FMIX{{{','.join(weight_statements)}}};"

    # write to file
    with open('test_nex.nex', "w") as nex_file:
        nex_file.write("#nexus\nbegin models;\n")

        for line in out_freqs:
            nex_file.write(" ".join(line) + "\n")

        nex_file.write(weight_line + "\n")

        nex_file.write("end;")

    return 'test_nex.nex'


def run_iqtrees(trees, alignment_address, avg_alpha, model, nexus_file, cores, leaves):
    total_tree_count = len(trees)

    i = 1
    for tree in trees:

        # render image of stitched-together-tree
        tree.render(f"test_{i}.png")
        tree.write(format=1, outfile=f"test_{i}.tree")

        iqtree_cmd = [
            "iqtree2",
            "-s", alignment_address,
            "--tree-fix", f"test_{i}.tree",
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

        print(f"Running iqtree funDi on Tree {i} out of {total_tree_count}...")
        subprocess.run(iqtree_cmd)
        i += 1


def generate_summary(tree_count):

    # get tree_properties
    tree_properties = {}

    # TODO: 
    # loop like
    # for i in range(tree_count):
    # 1-based-index so no
    i = 1
    while i <= tree_count:

        # TODO: read up on tuples
        # a list of three tree_properties, fundi log-likelihood, rho value and central branch length
        attribute = [0, 0, 0]

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

        tree_properties[i] = attribute
        i += 1

    # get the best tree based on funDi log-likelihood
    # TODO: this could be simpler
    best_tree_index = max(tree_properties, key=lambda key: tree_properties[key][0])
    best_tree = None

    # TODO: make use of the .treefile instead of the .iqtree file
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

    # we're overwriting the initial .png image,
    # but only for the best tree

    # style the tree
    tree_style = TreeStyle()
    tree_style.show_leaf_name = True
    tree_style.show_branch_length = True
    tree_style.branch_vertical_margin = 10

    # give each internal node an explicit name
    for node in best_tree.traverse():

        if not node.is_leaf():
            # TODO: try empty strings
            node.name = "node"

    best_tree.render(file_name=f"test_{best_tree_index}.png", tree_style=tree_style, units="px", w=800, h=1000)
    # TODO: root the best_tree, reference tree outgroup

    # print to summary file
    print("Generating summary...")

    with open("test_summary.txt", "w") as summary_file:

        summary_file.write(
            f"Tree {best_tree_index} has the largest funDi log-likelihood of {tree_properties[best_tree_index][0]}.\n"
            f"rho: {tree_properties[best_tree_index][1]}.\n"
            f"Central branch length: {tree_properties[best_tree_index][2]}\n"
        )

        # print tree with branch lengths
        # TODO: look into documentation of .get_ascii function
        summary_file.write(f"{best_tree.get_ascii(attributes=['name', 'dist'], show_internal=True)}\n\n")
        summary_file.write(f"See \"test_{best_tree_index}.png\" for a tree illustration.\n\n")

        # TODO: sort trees by log likelihood
        for tree, attribute in tree_properties.items():
            summary_file.write(f"funDi Log-likelihood of the tree {tree}: {attribute[0]}; "
                               f"rho: {attribute[1]}; "
                               f"central branch length: {attribute[2]}\n")

    print("Summary generated under 'test_summary.txt'.")


def main(args):
    try:
        master_tree = Tree(args.tree)
        # TODO: remove "redundant" variable names for arg arguments
        #  NOTE: I would have to go against this again, implementing this caused a lot of error because you'd have to
        #        change entry by entry. This potentially could mean a lot of headache if we change anything with args.
        alignment_address = args.alignment
        alignment = validate_alignment(master_tree, alignment_address)
        defined_groups = validate_def_file(master_tree, args.definition)
        model = args.mixture_model
        nexus_address = args.nexus
        cores = args.cores
        a_tree, b_tree = get_ref_subtrees(master_tree, defined_groups)

        # write subtrees into newick files
        a_tree.write(format=1, outfile="test_a.tree")
        b_tree.write(format=1, outfile="test_b.tree")

        # split alignment into two sub-alignments
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

        # check if the new trees generated by iqtree have the same topology
        a_tree = conform_iqtree_tree("test_a.treefile", a_tree)
        b_tree = conform_iqtree_tree("test_b.treefile", b_tree)

        trees = []
        denominator = int(1 / float(args.increment))
        proportions = [x / denominator for x in range(1, denominator)]
        a_branch = a_tree.get_children()[0].dist + a_tree.get_children()[1].dist
        b_branch = b_tree.get_children()[0].dist + b_tree.get_children()[1].dist

        print(f"branch a: {a_branch}, branch b: {b_branch}\n")

        # get alpha and beta from a cartesian product of proportions
        for alpha, beta in itertools.product(proportions, repeat=2):

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

        # NOTE: discuss purpose of custom nexus file with Hector
        if not nexus_address:
            avg_mixture_weights = calculate_weighted_average_mixture_weights(
                "test_a.iqtree", "test_b.iqtree", a_weight, b_weight
            )
            # generate nexus file:
            nexus_address = write_nexus_file(avg_mixture_weights, model)

        # second iqtree execution
        run_iqtrees(trees, alignment_address, avg_alpha, model, nexus_address, cores, a_tree.get_leaf_names())

        # To get the number of trees generated, we take number of proportions to the power of 2
        generate_summary(len(trees))

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
        "-i", "0.3"
    ]

    arguments = parser.parse_args()

    main(arguments)
