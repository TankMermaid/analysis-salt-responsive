#!/usr/bin/env python3
#
# author: scott olesen <swo@mit.edu>

# you should probably run tr ';_' '_-' on the taxa first

import argparse, sys
import dendropy as dp

def distance_sep(separator, ancestor, descendant):
    a_fields = ancestor.split(separator)
    d_fields = descendant.split(separator)

    if a_fields == d_fields:
        return 0
    elif len(a_fields) > len(d_fields):
        x = distance_sep(separator, descendant, ancestor)
        if x is None:
            return None
        else:
            return -x
    else:
        if d_fields[0: len(a_fields)] == a_fields:
            return len(d_fields) - len(a_fields)
        else:
            return None

def intermediates_sep(separator, ancestor, descendant):
    assert(distance_sep(separator, ancestor, descendant) > 0)
    a_fields = ancestor.split(separator)
    d_fields = descendant.split(separator)
    n_intermediates = len(d_fields) - len(a_fields) + 1
    return [separator.join(d_fields[0: len(a_fields) + i]) for i in range(1, n_intermediates)]

def nonleaf_singleton_children(node):
    return [c for c in node.child_nodes() if c.is_internal()]

if __name__ == '__main__':
    p = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('taxa', type=argparse.FileType('r'))
    p.add_argument('--separator', '-F', default='_')
    p.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'))
    p.add_argument('--suppress_inner', action='store_true', help='remove unifurcations?')
    args = p.parse_args()

    taxa = [l.rstrip() for l in args.taxa]

    distance = lambda ancestor, descendant: distance_sep(args.separator, ancestor, descendant)

    tree = dp.Tree()
    tree.seed_node.label = "Root"

    for taxon in taxa:
        # find the closest node on the tree
        closest = tree.seed_node
        closest_dist = distance(closest.label, taxon)
        for node in tree.levelorder_node_iter():
            dist = distance(node.label, taxon)
            if dist is not None and dist < closest_dist:
                closest = node
                closest_dist = dist

        # if this node is already on the tree, stop
        if closest_dist > 0:
            # figure out the necessary intermediates
            intermediate_labels = intermediates_sep(args.separator, closest.label, taxon)
            assert(taxon == intermediate_labels[-1])
            assert(closest.label not in intermediate_labels)
            for inter in intermediate_labels:
                closest = closest.add_child(dp.Node(edge_length=1, label=inter))

    # make sure everything comes to the maximum depth
    max_depth = max([leaf.distance_from_root() for leaf in tree.leaf_node_iter()])
    for leaf in tree.leaf_node_iter():
        leaf.edge.length += max_depth - leaf.distance_from_root()

    if args.suppress_inner:
        tree.suppress_unifurcations()

    tree.write_to_stream(args.output, 'newick', suppress_leaf_node_labels=False)
