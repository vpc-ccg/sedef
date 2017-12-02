# 786 

import pandas as pd
import numpy as np
import subprocess as sp
import sys, re, os, glob
from kswpython import kswalign, FastaReference
from collections import *

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

# a, b = "ATAGCTAGCTAGCAT", "AGCTAcCGCATCCC"
# aln = kswalign(a, b, 1, -2, 2, 1)
# print aln
# a, b, aln = process(a, b, aln)
# print a
# print aln
# print b

#--------

chrom1 = sys.argv[1]
chrom2 = sys.argv[2]
strand = sys.argv[3]
path = '{}_{}_{}.bed'.format(chrom1, chrom2, strand)
strand = '_' if strand == 'y' else '+'

frA = FastaReference("data/hg19/{}.fa".format(chrom1))
frB = FastaReference("data/hg19/{}.fa".format(chrom2))

def align(a, b):
    return kswalign(a, b, 5, -4, 40, 1)

# A-B is WGAC
def print_nicely(chromA, chromB, wgacA, wgacB, sedefA, sedefB):
    # fix rev comp later
    # seq_sedefA = frA.getSubSequence(chromA, sedefA[0], sedefA[1] - sedefA[0])
    # seq_sedefB = frB.getSubSequence(chromB, sedefB[0], sedefB[1] - sedefB[0])
    # print seq_sedefA
    # print seq_sedefB
    # aln_sedef = align(seq_sedefA, seq_sedefB)
    # print aln_sedef
    # seq_sedefA, seq_sedefB, aln_sedef = process(seq_sedefA, seq_sedefB, aln_sedef)

    seq_wgacA = frA.getSubSequence(chromA, wgacA[0], wgacA[1] - wgacA[0])
    seq_wgacB = frB.getSubSequence(chromB, wgacB[0], wgacB[1] - wgacB[0])
    aln_wgac = align(seq_wgacA, seq_wgacB)
    seq_wgacA, seq_wgacB, aln_wgac = process(seq_wgacA, seq_wgacB, aln_wgac)
    
    wgac_front = ''
    # sedefA_front, sedefB_front = '', ''
    # if sedefA[0] < wgacA[0]:
    #     wgac_front += '*' * (wgacA[0] - sedefA[0])
    # else:
    #     sedefA_front += '*' * (sedefA[0] - wgacA[0])
    # if sedefB[0] < wgacB[0]:
    #     if wgacB[0] - sedefB[0] > len(wgac_front):
    #         wgac_front += '*' * (wgacB[0] - sedefB[0] - len(wgac_front))
    #         sedefA_front += '*' * (wgacB[0] - sedefB[0] - len(wgac_front))
    #     else:
    #         sedefB_front += '*' * (len(wgac_front) - wgacB[0] - sedefB[0])
    # else:
    #     sedefB_front += '*' * (sedefB[0] - wgacB[0])

    # print 'SEDEF', sedefA_front + seq_sedefA
    blanks, smallA, smallB = 0, 0, 0

    for i in xrange(0, len(seq_wgacA), 100):
        print ' > ', seq_wgacA[i:i + 100]
        print ' + ', aln_wgac[i:i + 100]
        print ' > ', seq_wgacB[i:i + 100]
        blanks += sum(1 for c in aln_wgac[i:i + 100] if c == ' ')
        smallA += sum(1 for c in seq_wgacA[i:i + 100] if c.islower())
        smallB += sum(1 for c in seq_wgacB[i:i + 100] if c.islower())
        print '   ', i, blanks, smallA, smallB
    # print 'SEDEF', sedefB_front + seq_sedefB
    

print_nicely('chr22', 'chr22', 
    (32676644,  32788590),
    (44619626,  44731572), 0, 0)
exit(0)

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
system("bedtools pairtopair -a temp.bed -b <(cat {} | tr -d ,) -is -type both > temp_diff.bed".format(path))
print 'Paired WGAC and Sedef, size: {}'.format(system("wc -l temp_diff.bed"))
# exit(0)

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
    print_nicely(chrom1, chrom2, A[0:2], A[2:], B[0:2], B[2:])
    print
    #exit(0)

