# 786 

import pandas as pd
import numpy as np
import subprocess as sp
import sys, re, os, glob
from collections import *

def system(x):
    return sp.check_output(x, shell=True).strip()



# % normal bases
# % gaps per my calc
# % mis per my calc
# boundary Ms
# [via sedef] % completely covered cases

# for (auto &c: cigar) {
#     len += c.second;
#     if (prnaln) prnn("{}:{} ", c.first, c.second);
#     for (int i = 0; i < c.second; i++) {
#         assert(c.first != 'M' || ia < a.size());
#         assert(c.first != 'M' || ib < b.size());
#         if (c.first == 'M' && ceq(a[ia], b[ib])) {
#             get<2>(result) += "|";
#         } else {
#             get<2>(result) += " ";
#         }
#         if (c.first != 'D') get<1>(result) += b[ib++];
#         else                get<1>(result) += "-";
#         if (c.first != 'I') get<0>(result) += a[ia++];
#         else                get<0>(result) += "-";
#     }
# }

def perc(a, b):
    return 100.0 * a / float(b)

def process(a, b, cigar):
    cigar = [s for s in re.split('([MID;])', cigar) if s not in [';', '']]
    cigar = [(cigar[i + 1], int(cigar[i])) for i in xrange(0, len(cigar), 2)]
    #print cigar

    left = 0 if cigar[0][0] != 'M' else cigar[0][1]
    right = 0 if cigar[-1][0] != 'M' else cigar[-1][1]

    fuged = 0
    fuged_ext = 0
    gaps = 0
    mism = 0
    tlen = 0

    ia, ib = 0, 0
    for op, sz in cigar:
        tlen += sz
        for i in xrange(sz):
            if (ia < len(a) and a[ia].isupper()) or (ib < len(b) and b[ib].isupper()):
                fuged_ext += 1
            if op == 'M': 
                if a[ia].isupper() and b[ib].isupper(): 
                    fuged += 1
                if a[ia].upper() != b[ib].upper(): 
                    mism += 1
            else: 
                gaps += 1
            if op != 'D': ib += 1
            if op != 'I': ia += 1

    gaps = perc(gaps, tlen)
    mism = perc(mism, tlen)
    fpt = perc(fuged, tlen)

    return ((left, right), (gaps, mism), (fuged, fuged_ext, tlen))

y = []
n = 0
with open(sys.argv[1]) as f:
    next(f)
    for l in f:
        l = l.strip().split()
        q = 29
        n += 1
        y.append(process(l[q+1], l[q+2], l[q+0]) + (n,))

print len(y)

min_left  = min(y, key=lambda x: x[0][0])
min_right = min(y, key=lambda x: x[0][1])
print min_left[3], min_left[0][0]
print min_right[3], min_right[0][1]

max_gaps = max(y, key=lambda x: x[1][0])
max_mism = max(y, key=lambda x: x[1][1])
print max_gaps[3], max_gaps[1][0]
print max_mism[3], max_mism[1][1]

min_fug = min(y, key=lambda x: x[2][0])
min_fuge = min(y, key=lambda x: x[2][1])
print min_fug[3], min_fug[2][0]
print min_fuge[3], min_fug[2][1]
print sum(1 for s in y if s[2][1] < s[2][2]/10)

exit(0)



#--------

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

