from unittest import TestCase
from Mapper import fix_topology, get_info, get_averages, write_nexus_file
from ete3 import Tree


class Test(TestCase):
    master_tree = Tree("../data/Hector/master.tree")
    a_tree = Tree("../data/Hector/a.tree")
    b_tree = Tree("../data/Hector/b.tree")
    sub_a_tree = master_tree.get_children()[0]
    sub_b_tree = master_tree.get_children()[1]
    iqtree_file = "../data/Dandan/toy.subset1.aln.iqtree"

    def test_fix_topology_unrooted_scenario_1(self):
        new_tree = fix_topology(self.a_tree, self.sub_a_tree)
        self.assertEqual(0, new_tree.robinson_foulds(self.sub_a_tree)[0], "Topologies don't match!")

    def test_fix_topology_unrooted_scenario_2(self):
        alt_b_tree = Tree("(((7:0.705912,8:0.875553)1:0.313682,(9:0.227179,10:0.22797)1:0.403535)1:0.140418,6:0.498266,"
                          "(4:0.268157,5:0.285967)1:0.231266);")
        new_tree = fix_topology(alt_b_tree, self.sub_b_tree)
        self.assertEqual(0, new_tree.robinson_foulds(self.sub_b_tree)[0], "Topologies don't match!")

    def test_fix_topology_unrooted_scenario_3(self):
        new_tree = fix_topology(self.b_tree, self.sub_b_tree)
        self.assertEqual(0, new_tree.robinson_foulds(self.sub_b_tree)[0], "Topologies don't match!")

    def test_fix_topology_rooted_scenario_1(self):
        alt_a_tree = self.a_tree.copy("deepcopy")
        alt_a_tree.set_outgroup("1")
        new_tree = fix_topology(alt_a_tree, self.sub_a_tree)
        self.assertEqual(0, new_tree.robinson_foulds(self.sub_a_tree)[0], "Topologies don't match!")

    def test_fix_topology_rooted_scenario_2(self):
        alt_a_tree = self.a_tree.copy("deepcopy")
        alt_a_tree.set_outgroup("1")
        alt_a_tree.set_outgroup(alt_a_tree.get_common_ancestor("3", "23"))
        new_tree = fix_topology(alt_a_tree, self.sub_a_tree)
        self.assertEqual(0, new_tree.robinson_foulds(self.sub_a_tree)[0], "Topologies don't match!")

    def test_fix_topology_rooted_scenario_3(self):
        alt_b_tree = self.b_tree.copy("deepcopy")
        alt_b_tree.set_outgroup("4")
        new_tree = fix_topology(alt_b_tree, self.sub_b_tree)
        self.assertEqual(0, new_tree.robinson_foulds(self.sub_a_tree)[0], "Topologies don't match!")

    def test_get_info(self):
        ref_weights = {"1": 0.0371, "2": 0.3617, "3": 0.0136, "4": 0.1076, "5": 0.1450, "6": 0.0902, "7": 0.1195,
                       "8": 0.0693, "9": 0.0560, "10": 0.0000}
        ref_alpha = 0.5975

        weights, alpha = get_info(self.iqtree_file)

        self.assertEqual(ref_weights, weights, "Weights don't match!")
        self.assertEqual(ref_alpha, alpha, "Alphas don't match!")
