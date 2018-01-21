# 786 

import pandas as pd
import numpy as np
import subprocess as sp
import sys, re, os, glob
from collections import *
#from kswpython import kswalign, FastaReference

def system(x):
    return sp.check_output(x, shell=True, executable="/bin/bash").strip()

def process(a, b, cigar):
    cigar = [s for s in re.split('([MID;])', cigar) if s not in [';', '']]
    cigar = [(cigar[i + 1], int(cigar[i])) for i in xrange(0, len(cigar), 2)]
    sa, sb, gp = '', '', ''
    for op, sz in cigar:
        if op == 'I':  
            sa += '-' * sz
            gp += ' ' * sz
        elif op == 'D':   
            sb += '-' * sz
            gp += ' ' * sz    
        else:
            gp += ''.join(' |'[a[i].lower() == b[i].lower()] for i in xrange(sz))
        if op != 'D': 
            sb += b[:sz]
            b = b[sz:]
        if op != 'I':
            sa += a[:sz]
            a = a[sz:]
    return sa, sb, gp

#----------------------------------------------------------------------------------------

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
print ':: Loaded {} hits from WGAC'.format(df.shape[0])

hits = {}
with open('temp.bed', 'w') as f:
    for _, r in df.iterrows():
        print >>f, '\t'.join(map(str, [
            r.chrom, r.chromStart, r.chromEnd,
            r.otherChrom, r.otherStart, r.otherEnd,
            r.alignfile, 0, '+', '-' if strand == '_' else '+'
        ]))
        hits[(r.chromStart, r.chromEnd, r.otherStart, r.otherEnd)] = list()

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

system("bedtools pairtopair -a temp.bed -b <(cat {} | tr -d ,) -is -type both > temp_diff.bed".format(path))
print ':: After bedtools we have {} hits to process'.format(system("wc -l temp_diff.bed"))

with open('temp_diff.bed') as f:
    for l in f:
        l = l.strip().split()
        q = 0
        A = (int(l[q+1]), int(l[q+2]), int(l[q+4]), int(l[q+5])) # WGAC
        q = 10
        B = (int(l[q+1]), int(l[q+2]), int(l[q+4]), int(l[q+5])) # SEDEF
        hits[A].append((diff(A, B), A, B))
        if (A[2],A[3],A[0],A[1]) in hits:
            hits[(A[2],A[3],A[0],A[1])].append((diff(A, B), A, B, l[6]))

tm = sum(1 for k, vs in hits.items() if len(vs) == 0)
print ':: Missed {} hits ({:.1f}%)'.format(tm, 100.0*tm/df.shape[0])
for k, vs in sorted(hits.items()):
    if len(vs) == 0:
        print '   -- missed {}: {}'.format(vs[3], k)
haha = []
for k, vs in hits.items():
    if len(vs) == 0: continue
    ok = any(v[0][0] == 100 and v[0][2] == 100 for v in vs)
    if not ok: 
        haha += [max(vs, key=lambda x: x[0][0] + x[0][2])]

tm = len(haha)
print ':: Partial {} hits ({:.1f}%)'.format(tm, 100.0*tm/df.shape[0])
for ((p1, n1, p2, n2), A, B) in sorted(haha):
    print '   -- partial: {:.1f}% and {:.1f}% ({} and {}) -- {} to {}'.format(p1, p2, n1, n2, A, B)
    # print_nicely(chrom1, chrom2, A[0:2], A[2:], B[0:2], B[2:])
    # print
    #exit(0)

