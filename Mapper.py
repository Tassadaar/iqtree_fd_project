import sys
import copy
import argparse
from itertools import product
import re
from ete3 import Tree, TreeNode


def fix_topology(input_tree, reference_tree):
    tree_copy = copy.deepcopy(input_tree)
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


def get_info(in_file):
    weights = None
    alpha = None

    with open(in_file, "r") as file:

        for line in file:

            if "Mixture weights:" in line:
                words = line.split()
                weights = words[2:]
                print(weights)

            if "Gamma shape alpha:" in line:
                words = line.split()
                alpha = float(words[3])
                print(alpha)
                break

    return weights, alpha


def get_averages():
    a_log = "data/Dandan/toy.subset1.aln.log"
    b_log = "data/Dandan/toy.subset2.aln.log"

    a_weights, a_alpha = get_info(a_log)
    b_weights, b_alpha = get_info(b_log)

    avg_weights = [(float(a_weight) + float(b_weight)) / 2 for a_weight, b_weight in zip(a_weights, b_weights)]
    avg_alpha = (a_alpha + b_alpha) / 2

    return avg_weights, avg_alpha


def main(args):

    try:
        master_tree = Tree(args.master_tree)

        assert len(master_tree.get_children()) == 2, f"Master tree must be rooted!\n {master_tree}"

        denominator = int(1 / args.increment)
        a_tree = Tree(args.a_tree)
        b_tree = Tree(args.b_tree)

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
            new_master_tree.add_child(copy.deepcopy(new_a_tree))
            new_master_tree.add_child(copy.deepcopy(new_b_tree))

            assert new_master_tree.robinson_foulds(master_tree)[0] == 0, "The new master tree is not the same!"

            trees.append(new_master_tree)

        # print number of trees generated
        print(f"{len(trees)} tree were generated")

    except AssertionError as e:
        print(f"Oops! {e}")

    except NameError as e:
        print(f"Panda-monium! {e}")

    finally:
        print("\nSystem exiting...")
        sys.exit()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Tree mapper")

    parser.add_argument("-m", "--master_tree", required=True)
    parser.add_argument("-a", "--a_tree", required=True)
    parser.add_argument("-b", "--b_tree", required=True)
    parser.add_argument("-i", "--increment", required=False, default=0.1)

    # emulating commandline arguments for development
    sys.argv = [
        "Mapper.py",
        "-m", "data/Dandan/rooted_toy.newick",
        "-a", "data/Dandan/toy.subset1.newick",
        "-b", "data/Dandan/toy.subset2.newick"
    ]

    arguments = parser.parse_args()
    main(arguments)
