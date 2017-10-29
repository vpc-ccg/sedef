# 786 

import pandas as pd
import numpy as np
import subprocess as sp
import sys, re, os, glob
from collections import *

def system(x):
    return sp.check_output(x, shell=True).strip()

chrom1 = sys.argv[1]
chrom2 = sys.argv[2]
strand = sys.argv[3]
path = '{}_{}_{}.bed'.format(chrom1, chrom2, strand)
strand = '_' if strand == 'y' else '+'

df = pd.read_table("data/GRCh37GenomicSuperDup.tab")
if chrom1 != chrom2 or strand == '_':
    df = df[(df.chrom == chrom1) & (df.otherChrom == chrom2) & (df.strand == strand)]
else:
    df = df[(df.chrom == chrom1) & (df.otherChrom == chrom2) & (df.strand == strand) & (df.chromStart < df.otherStart)]
df['chromSize'] = df.chromEnd - df.chromStart
print 'Loaded {} from WGAC'.format(df.shape[0])

hits = {}
with open('temp.bed', 'w') as f:
    for _, r in df.iterrows():
        print >>f, '\t'.join(map(str, [
            r.chrom, r.chromStart, r.chromEnd,
            r.otherChrom, r.otherStart, r.otherEnd
        ]))
        hits[(r.chromStart, r.chromEnd, r.otherStart, r.otherEnd)] = list()
print 'Wrote {} from WGAC to temp.bed'.format(df.shape[0])

def diff(X, Y):
    sA, eA = X[:2]
    sB, eB = Y[:2]
    
    x = max(sA, sB)
    y = min(eA, eB)
    d = max(y - x, 0) 
    wA = d
    dA = 100.0 * d / float(eA - sA)

    sA, eA = X[2:4]
    sB, eB = Y[2:4]
    
    x = max(sA, sB)
    y = min(eA, eB)
    d = max(y - x, 0)
    wB = d
    dB = 100.0 * d / float(eA - sA)

    return (dA, wA, dB, wB)

# bedtools
system("bedtools pairtopair -a temp.bed -b {} -is -type both > temp_diff.bed".format(path))
print 'Paired WGAC and Sedef, size: {}'.format(system("wc -l temp_diff.bed"))

with open('temp_diff.bed') as f:
    for l in f:
        l = l.strip().split()
        A = (int(l[1]), int(l[2]), int(l[4]), int(l[5]))
        B = (int(l[7]), int(l[8]), int(l[10]), int(l[11]))
        hits[A].append((diff(A, B), A, B))
        if (A[2],A[3],A[0],A[1]) in hits:
            hits[(A[2],A[3],A[0],A[1])].append((diff(A, B), A, B))

tm = sum(1 for k, vs in hits.items() if len(vs) == 0)
print 'total missed: {} ({:.1f}%)'.format(tm,100.0*tm/df.shape[0])
for k, vs in sorted(hits.items()):
    if len(vs) == 0:
        print '  missed', k
haha = []
for k, vs in hits.items():
    if len(vs) == 0: continue
    ok = any(v[0][0] == 100 and v[0][2] == 100 for v in vs)
    if not ok: 
        haha += [max(vs, key=lambda x: x[0][0] + x[0][2])]
for ((p1, n1, p2, n2), A, B) in sorted(haha):
    print 'partial: {:.1f}% and {:.1f}% ({} and {}) -- {} to {}'.format(p1, p2, n1, n2, A, B)

