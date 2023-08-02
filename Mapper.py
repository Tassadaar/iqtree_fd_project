import copy
import sys
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


def main():

    try:
        master_tree = Tree("((1:0.496516432,((21:0.2,(23:0.4,24:0.22):0.3):0.2750396589,3:0.4785155756):0.14105667032):0.1,"
                       "(((4:0.2864272562,5:0.2870235802):0.2266968989,6:0.5159850109):0.1365803773,"
                       "((7:0.7543316696,8:0.9338300134):0.34184149945,(9:0.2355059988,10:0.2518571253)"
                       ":0.3972088265):0.4108541642):0.01);")
        assert len(master_tree.get_children()) == 2, "Master tree must be rooted!"
    except AssertionError as e:
        print(f"Oops! {e}")
        sys.exit()

    a_tree = Tree("(1:0.6121762192,(21:0.1881984527,(23:0.4160042306,24:0.2186291776):0.3451463477)"
                      ":0.2916847347,3:0.5092916847);")
    b_tree = Tree("(4:0.2681573187,5:0.2859668232,(6:0.4982657718,((7:0.7059120269,8:0.8755530768):0.3136818250,"
                      "(9:0.2271793309,10:0.2279696348):0.4035351373):0.5616725078):0.2312659556);")
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
    main()
