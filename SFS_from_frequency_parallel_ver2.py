#!/usr/bin/env python3
"""
open the frequency file read backward in time until n individuals all coalesce.
As we go, we will record the number of leaves for each node
because the mutations happen uniformly in time, histogram of the number
of leaves of each branch will just be basically like SFS."""
import numpy as np
from multiprocessing import Pool
import sys

L = int(sys.argv[1])
N = int(sys.argv[2])
s = float(sys.argv[3])
m = float(sys.argv[4])
tfinal = int(sys.argv[5])
nbase = int(sys.argv[6])
N_SFS = int(sys.argv[7]) # number of coalescent simulation we run.

freq_file = open('L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}.txt'.format(L, N, s, m, tfinal))
lines = freq_file.readlines()

def get_parent(Ne_cumsum, Ne_parent, i):
    Ne_cumsum_parent = np.cumsum(Ne_parent)
    j = np.searchsorted(Ne_cumsum, i) + 1 # idx of first Ne_cumsum > i
    try:
        if j == 0 or Ne_cumsum_parent[j - 2] == Ne_cumsum_parent[j + 1]:
            return np.random.choice(range(Ne_cumsum_parent[1])
                , p = np.append((1 - m / 2) * np.ones(Ne_parent[0])
                , m / 2 * np.ones(Ne_parent[1]))
                / ((1 - m / 2) * Ne_parent[0] + m / 2 * Ne_parent[1]))
        elif j == 1:
            return np.random.choice(np.arange(0, Ne_cumsum_parent[j + 1])
                , p = np.concatenate((m / 2 * np.ones(Ne_parent[j - 1])
                , (1 - m) * np.ones(Ne_parent[j])
                , m / 2 * np.ones(Ne_parent[j + 1])))
                / (m / 2 * Ne_parent[j - 1]
                + (1 - m) * Ne_parent[j]
                + m / 2 * Ne_parent[j + 1]))
        else:
            return np.random.choice(np.arange(Ne_cumsum_parent[j - 2], Ne_cumsum_parent[j + 1])
                , p = np.concatenate((m / 2 * np.ones(Ne_parent[j - 1])
                , (1 - m) * np.ones(Ne_parent[j])
                , m / 2 * np.ones(Ne_parent[j + 1])))
                / (m / 2 * Ne_parent[j - 1]
                + (1 - m) * Ne_parent[j]
                + m / 2 * Ne_parent[j + 1]))
    except ValueError:
        print(j, Ne_cumsum_parent[1], Ne_cumsum_parent[j + 1],
              Ne_cumsum_parent[j - 2])
        raise

def runner(idx):
    Ne_0 = [int(N) for N in lines[-1].split()]
    Ne_cumsum_0 = np.cumsum(Ne_0)
    n = nbase
    SFS = np.zeros(n)

    Ne = Ne_0
    Ne_cumsum = Ne_cumsum_0
    if n < sum(Ne):
        individuals = np.random.choice(range(sum(Ne)), size = n, replace = False)
    else:
        individuals = range(sum(Ne))
        n = sum(Ne)
    unique, leaf_counts = np.unique(individuals, return_counts = True)
    hist, bin_edges = np.histogram(leaf_counts, bins = np.arange(1, n + 2))

    line_num = -1
    while (len(individuals) > 1) and (line_num > -len(lines)):
        line_num -= 1
        Ne_parent = [int(N) for N in lines[line_num].split()]
        individuals2 = []
        for i in individuals:
            individuals2.append(get_parent(Ne_cumsum, Ne_parent, i))
        individuals2 = np.repeat(individuals2, leaf_counts)
        unique, leaf_counts = np.unique(individuals2, return_counts = True)
        hist, bin_edges = np.histogram(leaf_counts, bins = np.arange(1, n + 2))
        SFS += hist
        individuals = unique
        print('Run #', idx, len(individuals2), len(unique))
    #    print(individuals)
        Ne = Ne_parent
        Ne_cumsum = np.cumsum(Ne_parent)
    if np.mod(idx, 100) == 0:
        np.savetxt('SFS_L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_nsample={}_{}.txt'.format(L,
           N, s, m, tfinal, n, idx), SFS)
    return SFS

if __name__ == '__main__':

        # this is the true sample number in case Ne < nbase.
    #print(individuals)
    p = Pool(10)
    SFS_items = p.map(runner, range(N_SFS))
    SFS = np.sum(SFS_items, axis=0)

    SFS /= N_SFS
    np.savetxt('SFS_L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_nsample={}_navg={}.txt'.format(L,
               N, s, m, tfinal, nbase, N_SFS), SFS)

