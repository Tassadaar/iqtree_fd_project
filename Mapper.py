import sys
import copy
import argparse
from ete3 import Tree, TreeNode


def fix_topology(input_tree, result_tree):
    new_tree = copy.deepcopy(input_tree)
    new_leaves = [child.get_leaf_names() for child in new_tree.get_children()]
    result_leaves = [child.get_leaf_names() for child in result_tree.get_children()]

    for leaf_group in new_leaves:

        if leaf_group not in result_leaves:
            continue

        # scenario 1: leaf group is valid but only one taxon
        if len(leaf_group) == 1:
            new_tree.set_outgroup(leaf_group.pop())
            return new_tree

        # scenario 2: leaf group is valid but multiple taxa
        leaf_nodes = [new_tree.get_leaves_by_name(leaf_name) for leaf_name in leaf_group]
        outgroup = leaf_nodes.pop().get_common_ancestor(leaf_nodes)
        new_tree.set_outgroup(outgroup)
        return new_tree

    # scenario 3: none of the leaf groups of root is valid
    def find_outgroup(node):
        parent = node.get_ancestors()[0]

        if parent.get_leaf_names() in result_leaves:
            return parent
        else:
            return find_outgroup(parent)

    new_tree.set_outgroup(find_outgroup(new_tree.get_farthest_leaf()[0]))
    return new_tree


def main(args):

    try:
        master_tree = Tree(args.master_tree)
        assert len(master_tree.get_children()) == 2, "Master tree must be rooted!"
    except AssertionError as e:
        print(f"Oops! {e}")
        sys.exit()

    a_tree = Tree(args.a_tree)
    b_tree = Tree(args.b_tree)
    sub_a_tree = master_tree.get_children()[0]
    sub_b_tree = master_tree.get_children()[1]

    new_a_tree = fix_topology(a_tree, sub_a_tree)
    new_b_tree = fix_topology(b_tree, sub_b_tree)
    new_master_tree = TreeNode()
    new_master_tree.add_child(new_a_tree)
    new_master_tree.add_child(new_b_tree)

    print(new_master_tree.robinson_foulds(master_tree)[0])
    print(new_master_tree)


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
