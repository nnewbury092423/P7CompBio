#!/usr/bin/env python3
ERROR_GTR_PARAMS_FILE = "Invalid GTR parameters file"
ERROR_MSA = "Invalid multiple sequence alignment"
ERROR_SCORE = "Invalid log-likelihood score. Must be a negative float"

def likelihood(tree, seqs, gtr_probs, gtr_rates):
    '''
    This function uses Felsenstein's tree-pruning algorithm to compute the log-likelihood score of a tree given an MSA
    :param tree: A TreeSwift Tree object representing the phylogenetic tree
    :param seqs: A dictionary where keys are labels corresponding to the labels of the leaves in ``tree`` and values are sequences (strings)
    :param gtr_probs: The GTR stationary probabilities as a list [prob_A, prob_C, prob_G, prob_T]
    :param gtr_rates: The GTR transition rates as a list [rate_CT, rate_AT, rate_GT, rate_AC, rate_CG, rate_AG]
    :return: The log-likelihood score
    '''
    prob_A, prob_C, prob_G, prob_T = gtr_probs # You can use these if it's more convenient
    rate_CT, rate_AT, rate_GT, rate_AC, rate_CG, rate_AG = gtr_rates # You can use these if it's more convenient
    log_likelihood_score = 0. # You can modify this variable, or you can use your own; either way is fine
    # TODO Your code here
    return log_likelihood_score

def read_FASTA(filename):
    '''
    This function reads a FASTA file from a file and returns a dictionary mapping identifiers to sequences
    '''
    stream = open(filename); seqs = {}; name = None; seq = ''
    for line in stream:
        l = line.strip()
        if len(l) == 0:
            continue
        if l[0] == '>':
            if name is not None:
                assert len(seq) != 0, "Malformed FASTA"
                seqs[name] = seq
            name = l[1:]
            assert name not in seqs, "Duplicate sequence ID: %s" % name
            seq = ''
        else:
            seq += l
    assert name is not None and len(seq) != 0, "Malformed FASTA"
    seqs[name] = seq; stream.close()
    return seqs

if __name__ == "__main__":
    '''
    This is how we handle loading the input dataset, running your functions, and outputting the results
    '''
    # parse user arguments
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tree', required=False, type=str, default='stdin', help="Input Tree File (Newick)")
    parser.add_argument('-p', '--gtr_params', required=True, type=str, help="GTR Parameters File")
    parser.add_argument('-s', '--seqs', required=True, type=str, help="Multiple Sequence Alignment (FASTA)")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output Log-Likelihood Score")
    args = parser.parse_args()

    # load input tree
    from treeswift import read_tree_newick
    if args.tree == 'stdin':
        from sys import stdin
        tree = read_tree_newick(stdin)
    else:
        tree = read_tree_newick(args.tree)

    # load GTR parameters
    try:
        gtr_probs, gtr_rates = [[float(e.strip()) for e in l.strip().split()] for l in open(args.gtr_params)]
    except:
        raise ValueError(ERROR_GTR_PARAMS_FILE)

    # load MSA
    try:
        seqs = read_FASTA(args.seqs)
    except:
        raise ValueError(ERROR_MSA)
    for l in tree.traverse_leaves():
        assert l.label in seqs, "Missing sequence for: %s" % l.label

    # run student code and output
    score = likelihood(tree, seqs, gtr_probs, gtr_rates)
    if (not isinstance(score,float) and not isinstance(score,int)) or score > 0:
        raise ValueError(ERROR_SCORE)
    if args.output == 'stdout':
        from sys import stdout as outfile
    else:
        outfile = open(args.output,'w')
    outfile.write('%s\n' % str(score)); outfile.close()
