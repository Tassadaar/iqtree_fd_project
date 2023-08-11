from unittest import TestCase
from Mapper import fix_topology
from ete3 import Tree


class Test(TestCase):
    master_tree = Tree("((1:0.496516432,((21:0.2,(23:0.4,24:0.22):0.3):0.2750396589,3:0.4785155756):0.14105667032):0.1,"
                       "(((4:0.2864272562,5:0.2870235802):0.2266968989,6:0.5159850109):0.1365803773,"
                       "((7:0.7543316696,8:0.9338300134):0.34184149945,(9:0.2355059988,10:0.2518571253)"
                       ":0.3972088265):0.4108541642):0.01);")
    a_tree = Tree("(1:0.6121762192,(21:0.1881984527,(23:0.4160042306,24:0.2186291776):0.3451463477)"
                  ":0.2916847347,3:0.5092916847);")
    b_tree = Tree("(4:0.2681573187,5:0.2859668232,(6:0.4982657718,((7:0.7059120269,8:0.8755530768):0.3136818250,"
                  "(9:0.2271793309,10:0.2279696348):0.4035351373):0.5616725078):0.2312659556);")
    sub_a_tree = master_tree.get_children()[0]
    sub_b_tree = master_tree.get_children()[1]

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

    def test_fix_topology_rooted_single_outgroup(self):
        pass

    def test_fix_topology_rooted_multiple_outgroup(self):
        pass
