""" This script transforms an anonymized tree into a tree with random uniprot ids that overlap in taxid. """

import sys
import re
import random

tree_file = sys.argv[1]
with open(tree_file) as f:
    tree_contents = f.read()

n_leaves = len(re.findall("L\d{5}", tree_contents))

tree_contents = re.sub("L\d{5}", "{}", tree_contents)

uniprots = []
skip = random.randint(0, 100)
with open("../taxids/taxids_small.txt") as f:
    for l, line in enumerate(f):
        if l < skip:
            continue
        if l - skip > n_leaves:
            break
        uniprot, _, taxid = line.strip().split("\t")
        uniprots.append(uniprot)

new_tree = tree_contents.format(*uniprots)
with open(tree_file.replace("anon", "random"), "w+") as f:
    f.write(new_tree)
