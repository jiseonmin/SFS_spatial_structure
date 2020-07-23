#!/usr/bin/env python3

# Debug the factor difference. (neutral theory : 2NU/f)
import numpy as np
import matplotlib.pyplot as plt
L = 2000
N = 100
s = 0.05
m = 0.25
tfinal = 10000
n = 10000
nSFS = 1000

# Plot for different s
f = np.arange(1, n + 1) / n
plt.figure(figsize = (8, 6))

for m in [0.25, 0.15, 0.1]:
    SFS = np.loadtxt('SFS_L={}_N={}_s={:.6f}_m={:.6f}_tfinal={}_nsample={}_sample_middle_avg.txt'.format(L, 
               N, s, m, tfinal, n))

    plt.loglog(f, SFS / nSFS, label = 'm = {:.2f}'.format(m), linewidth = 0.9)


plt.loglog(f[:1000], 1 / f[:1000] 
* f[0] * SFS[0] / nSFS, label = '$c/f$', linewidth = 2)
##fmiddle = f[50:1000]
plt.loglog(f[:1000],
           1 / f[:1000] ** 2 * f[0] ** 2 *
           SFS[0] / nSFS, label = '$c_2/f^2$', linewidth = 2)
plt.xlabel('f')
plt.ylabel('n(f)')
plt.legend()
plt.title('N = {}, t = {}, s = {:.2f}, sample middle'.format(N, tfinal, s))