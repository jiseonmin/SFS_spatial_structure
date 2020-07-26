#!/usr/bin/env python3
"""
Coalescent simulation of individuals in L demes (N per deme).
As we go, we will record the number of leaves for each node. Then we sample 
some of the leaf counts with rate Un, since a neutral mutation arises as a 
Poisson process. We append the counts until we reach a single common ancestor.
The histogram of the counts is SFS """
import numpy as np
from multiprocessing import Pool
import sys

L = int(sys.argv[1])
N = int(sys.argv[2])
m = float(sys.argv[3])
tfinal = int(sys.argv[4])
Un = float(sys.argv[5])
nbase = int(sys.argv[6])
N_SFS = int(sys.argv[7]) # number of coalescent simulation we run.

def get_parent(Ne_cumsum, Ne, i):
#    Ne_cumsum_parent = np.cumsum(Ne_parent)
    j = np.searchsorted(Ne_cumsum, i) + 1 # idx of first Ne_cumsum > i
    try:
        if j == 0 or Ne_cumsum[j - 2] == Ne_cumsum[j + 1]:
            return np.random.choice(range(Ne_cumsum[1])
                , p = np.append((1 - m / 2) * np.ones(Ne[0])
                , m / 2 * np.ones(Ne[1]))
                / ((1 - m / 2) * Ne[0] + m / 2 * Ne[1]))
        elif j == 1:
            return np.random.choice(np.arange(0, Ne_cumsum[j + 1])
                , p = np.concatenate((m / 2 * np.ones(Ne[j - 1])
                , (1 - m) * np.ones(Ne[j])
                , m / 2 * np.ones(Ne[j + 1])))
                / (m / 2 * Ne[j - 1]
                + (1 - m) * Ne[j]
                + m / 2 * Ne[j + 1]))
        else:
            return np.random.choice(np.arange(Ne_cumsum[j - 2], Ne_cumsum[j + 1])
                , p = np.concatenate((m / 2 * np.ones(Ne[j - 1])
                , (1 - m) * np.ones(Ne[j])
                , m / 2 * np.ones(Ne[j + 1])))
                / (m / 2 * Ne[j - 1]
                + (1 - m) * Ne[j]
                + m / 2 * Ne[j + 1]))
    except ValueError:
        print(j, Ne_cumsum[1], Ne_cumsum[j + 1],
              Ne_cumsum[j - 2])
        raise

def runner(idx):
    Ne = L * np.ones(N)
    Ne_cumsum = np.cumsum(Ne)
    n = nbase
    SFS = np.zeros(n)

    if n < sum(Ne):
        individuals = np.random.choice(range(sum(Ne)), size = n, replace = False)
    else:
        individuals = range(sum(Ne))
        n = sum(Ne)
    unique, leaf_counts = np.unique(individuals, return_counts = True)
    hist, bin_edges = np.histogram(leaf_counts, bins = np.arange(1, n + 2))

    while (len(individuals) > 1):
#        Ne_parent = [int(N) for N in lines[line_num].split()]
        individuals2 = []
        for i in individuals:
            individuals2.append(get_parent(Ne_cumsum, Ne, i))
        individuals2 = np.repeat(individuals2, leaf_counts)
        unique, leaf_counts = np.unique(individuals2, return_counts = True)
        neutal_mut_counts = np.random.poisson(Un, len(unique))
        descendents_counts = np.repeat(leaf_counts, neutal_mut_counts)
        hist, bin_edges = np.histogram(descendents_counts, bins = np.arange(1, n + 2))
        SFS += hist
        individuals = unique
        print('Run #', idx, len(individuals2), len(unique))
    #    print(individuals)
    if np.mod(idx, 10) == 0:
        np.savetxt('SFS_neutral_L={}_N={}_m={:.6f}_tfinal={}_nsample={}_Un={:.6f}_{}.txt'.format(L,
           N, m, tfinal, n, Un, idx), SFS)
    return SFS

if __name__ == '__main__':

        # this is the true sample number in case Ne < nbase.
    #print(individuals)
    p = Pool(10)
    SFS_items = p.map(runner, range(N_SFS))
    SFS = np.sum(SFS_items, axis=0)

    SFS /= N_SFS
    np.savetxt('SFS_neutral_L={}_N={}_m={:.6f}_tfinal={}_nsample={}_Un={:.6f}_navg={}.txt'.format(L,
               N, m, tfinal, nbase, Un, N_SFS), SFS)

