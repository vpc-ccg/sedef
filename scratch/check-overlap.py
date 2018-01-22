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
path = sys.argv[4]
#path = '{}_{}_{}.bed'.format(chrom1, chrom2, strand)
strand = '_' if strand == 'y' else '+'

df = pd.read_table("data/GRCh37GenomicSuperDup.tab")
if chrom1 != chrom2 or strand == '_':
    df = df[(df.chrom == chrom1) & (df.otherChrom == chrom2) & (df.strand == strand)]
else:
    df = df[(df.chrom == chrom1) & (df.otherChrom == chrom2) & (df.strand == strand) & (df.chromStart < df.otherStart)]
df['chromSize'] = df.chromEnd - df.chromStart
print ':: Loaded {} hits from WGAC'.format(df.shape[0])

hits = {}
name_to_coor = {}
with open('temp.bed', 'w') as f:
    for _, r in df.iterrows():
        print >>f, '\t'.join(map(str, [
            r.chrom, r.chromStart, r.chromEnd,
            r.otherChrom, r.otherStart, r.otherEnd,
            r.alignfile, 0, '+', '-' if strand == '_' else '+'
        ]))
        hits[r.alignfile] = list()
        name_to_coor[r.alignfile] = (r.chromStart, r.chromEnd, r.otherStart, r.otherEnd)


def diff(wgac, sedef): # how much wgac is off sedef
    def overlap(sa, ea, sb, eb):
        return max(0, min(ea, eb) - max(sa, sb))

    def match(sW, eW, sS, eS):
        oo = overlap(sW, eW, sS, eS)
        wW = oo
        dW = 100.0 * wW / float(eW - sW)
        return (dW, wW, eW - sW)

    # left match
    sW, eW = wgac[:2]
    sS, eS = sedef[:2]

    r1  = match(*(wgac[:2] + sedef[:2]))
    r1 += match(*(wgac[2:4] + sedef[2:4]))

    r2  = match(*(wgac[2:4] + sedef[:2]))
    r2 += match(*(wgac[:2] + sedef[2:4]))

    return r1 if r1[0]+r1[3] > r2[0]+r2[3] else r2

system("bedtools pairtopair -a temp.bed -b <(cat {} | tr -d ,) -is -type both > temp_diff.bed".format(path))
print ':: After bedtools we have {} hits to process'.format(system("wc -l temp_diff.bed"))

with open('temp_diff.bed') as f:
    for l in f:
        l = l.strip().split()
        q = 0
        A = (int(l[q+1]), int(l[q+2]), int(l[q+4]), int(l[q+5])) # WGAC
        q = 10
        B = (int(l[q+1]), int(l[q+2]), int(l[q+4]), int(l[q+5])) # SEDEF
        name = l[6]
        hits[name].append((diff(A, B), A, B))

try:
    tm = sum(1 for k, vs in hits.items() if len(vs) == 0)
    print ':: Missed {} hits ({:.1f}%)'.format(tm, 100.0*tm/df.shape[0])
    for name, vs in sorted(hits.items()):
        if len(vs) == 0:
            print '   -- missed http://humanparalogy.gs.washington.edu/build37/{0}'.format(name)
            p = name_to_coor[name]
            p += (p[1] - p[0], p[3] - p[2])
            print '      wgac: {:11,}..{:11,} -> {:11,}..{:11,} ... len {:9,} -> {:9,} '.format(*p)

    haha = defaultdict(list)
    for k, vs in hits.items():
        if len(vs) == 0: continue
        ok = any(round(v[0][0], 2) >= 99.99 and round(v[0][3], 2) >= 99.99 for v in vs)
        if not ok: 
            haha[k] += vs # [max(vs, key=lambda x: x[0][0] + x[0][3])]

    tm = len(haha)
    print ':: Partial {} hits ({:.1f}%)'.format(tm, 100.0*tm/df.shape[0])
    for k in sorted(haha.keys(), key=lambda y: sum(yy[0][0] + yy[0][3] for yy in haha[y])):
        print '   -- partial http://humanparalogy.gs.washington.edu/build37/{0}'.format(k)
        for ((p1, n1, t1, p2, n2, t2), A, B) in sorted(haha[k]):
            print '      === {:.1f}% ({} of {}) and {:.1f}% ({} of {})'.format(p1, n1, t1, p2, n2, t2)
            A += (A[1] - A[0], A[3] - A[2])
            B += (B[1] - B[0], B[3] - B[2])
            print '          wgac:  {:11,}..{:11,} -> {:11,}..{:11,} ... len {:9,} -> {:9,} '.format(*A)
            print '          sedef: {:11,}..{:11,} -> {:11,}..{:11,} ... len {:9,} -> {:9,} '.format(*B)
except IOError:
    pass
