import sys
import copy
import argparse
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
        outgroup = leaf_group.pop().get_common_ancestor(leaf_group)
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


def main(args):

    try:
        master_tree = Tree(args.master_tree)

        assert len(master_tree.get_children()) == 2, "Master tree must be rooted!"

        a_tree = Tree(args.a_tree)
        b_tree = Tree(args.b_tree)

        for subtree in master_tree.get_children():
            if subtree.get_leaf_names() == a_tree.get_leaf_names():
                sub_a_tree = subtree
            if subtree.get_leaf_names() == b_tree.get_leaf_names():
                sub_b_tree = subtree

        new_a_tree = fix_topology(a_tree, sub_a_tree)
        new_b_tree = fix_topology(b_tree, sub_b_tree)

        # reconstruct master tree
        new_master_tree = TreeNode()
        new_master_tree.add_child(new_a_tree)
        new_master_tree.add_child(new_b_tree)

        assert new_master_tree.robinson_foulds(master_tree)[0] == 0, "The new master tree is not the same!"

        print(f"Reference master tree:{master_tree} \nSubtree A:{a_tree} \nSubtree B:{b_tree}"
              f"\nNew master tree:{new_master_tree}")

    except AssertionError as e:
        print(f"Oops! {e}")

    except NameError as e:
        print(f"Panda-monium! {e}")

    finally:
        sys.exit()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Tree mapper")

    parser.add_argument("-m", "--master_tree", required=True)
    parser.add_argument("-a", "--a_tree", required=True)
    parser.add_argument("-b", "--b_tree", required=True)

    # emulating commandline arguments for development
    sys.argv = [
        "Mapper.py",
        "-m", "master.tree",
        "-a", "a.tree",
        "-b", "b.tree"
    ]

    arguments = parser.parse_args()
    main(arguments)
