import sys
import dendropy
from dendropy.calculate import treecompare
sys.setrecursionlimit(10000)


def read_tree(fname, tns):
    return dendropy.Tree.get_from_path(fname, "newick", taxon_namespace=tns)

def main(tree_path_1, tree_path_2):
    tns = dendropy.TaxonNamespace()

    tree1 = read_tree(tree_path_1, tns)
    tree2 = read_tree(tree_path_2, tns)

    tree1.encode_bipartitions()
    tree2.encode_bipartitions()
    

    print("Number of leaves in tree 1:         ", len(tree1.leaf_nodes()))
    print("Number of leaves in tree 2:         ", len(tree2.leaf_nodes()))
    print("Unweighted Robinson-Fould distance: ", treecompare.symmetric_difference(tree1, tree2))
    print("Weighted Robinson-Fould distance:   ", treecompare.weighted_robinson_foulds_distance(tree1, tree2))
    print("Euclidean distance:                 ", treecompare.euclidean_distance(tree1, tree2))

if __name__ == "__main__":
    tree_path_1, tree_path_2 = sys.argv[1], sys.argv[2]

    main(tree_path_1, tree_path_2)
