import sys
import dendropy
from dendropy.calculate import treecompare

tree_path_1, tree_path_2 = sys.argv[1], sys.argv[2]

tns = dendropy.TaxonNamespace()
tree1 = dendropy.Tree.get_from_path(
        tree_path_1,
        "newick",
        taxon_namespace=tns)
tree2 = dendropy.Tree.get_from_path(
        tree_path_2,
        "newick",
        taxon_namespace=tns)

tree1.encode_bipartitions()
tree2.encode_bipartitions()

print("Unweighted Robinson-Fould distance: ", treecompare.symmetric_difference(tree1, tree2))
print("Weighted Robinson-Fould distance:   ", treecompare.weighted_robinson_foulds_distance(tree1, tree2))
print("Euclidean distance:                 ", treecompare.euclidean_distance(tree1, tree2))