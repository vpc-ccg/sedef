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

path = sys.argv[1]
if len(sys.argv) > 2:
    chrom1 = sys.argv[2]
    chrom2 = sys.argv[3]
    strand = sys.argv[4]
    strand = '_' if strand == 'y' else '+'
else:
    chrom1 = ''
#path = '{}_{}_{}.bed'.format(chrom1, chrom2, strand)

# sizes = {}
# with open('data/hg19/hg19.fa.fai') as f:
#     for l in f:
#         l = l.strip().split()
#         sizes[l[0]] = l[1]

df = pd.read_table("data/GRCh37GenomicSuperDup.tab")

if chrom1 != '':
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
        if '_' in r.chrom or '_' in r.otherChrom:
            continue
        print >>f, '\t'.join(map(str, [
            r.chrom, r.chromStart, r.chromEnd,
            r.otherChrom, r.otherStart, r.otherEnd,
            r.alignfile, 0, '+', '-' if r.strand == '_' else '+'
        ]))
        hits[r.alignfile] = list()
        name_to_coor[r.alignfile] = r #(r.chromStart, r.chromEnd, r.otherStart, r.otherEnd)


def diff(wgac, sedef): # how much wgac is off sedef
    def overlap(sa, ea, sb, eb):
        return max(0, min(ea, eb) - max(sa, sb))

    def match(sW, eW, sS, eS):
        oo = overlap(sW, eW, sS, eS)
        wW = oo
        dW = 100.0 * wW / float(eW - sW)
        return (dW, wW, eW - sW)

    r1  = match(*(wgac[1:3] + sedef[1:3]))
    r1 += match(*(wgac[4:6] + sedef[4:6]))

    r2  = match(*(wgac[4:6] + sedef[1:3]))
    r2 += match(*(wgac[1:3] + sedef[4:6]))

    return r1 if r1[0]+r1[3] > r2[0]+r2[3] else r2

system("bedtools pairtopair -a temp.bed -b <(cat {} | tr -d ,) -is -type both > temp_diff.bed".format(path))
print ':: After bedtools we have {} hits to process'.format(system("wc -l temp_diff.bed"))

with open('temp_diff.bed') as f:
    for l in f:
        l = l.strip().split()
        q = 0
        A = (l[q], int(l[q+1]), int(l[q+2]), l[q+3], int(l[q+4]), int(l[q+5])) # WGAC
        q = 10
        B = (l[q], int(l[q+1]), int(l[q+2]), l[q+3], int(l[q+4]), int(l[q+5])) # SEDEF
        name = l[6]
        hits[name].append((diff(A, B), A, B))

try:
    tm = sum(1 for k, vs in hits.items() if len(vs) == 0)
    print ':: Missed {} hits ({:.1f}%)'.format(tm, 100.0*tm/df.shape[0])
    for name, vs in sorted(hits.items()):
        if len(vs) == 0:
            print '   -- missed http://humanparalogy.gs.washington.edu/build37/{0}'.format(name)
            r = name_to_coor[name]
            print '      wgac: {:5} {:11,}..{:11,} -> {:5} {:11,}..{:11,} ... {} ... len {:9,} -> {:9,} '.format(
                r.chrom, r.chromStart, r.chromEnd, r.otherChrom, r.otherStart, r.otherEnd, r.strand,
                r.chromEnd - r.chromStart, r.otherEnd - r.otherStart
            )

    partials = defaultdict(list)
    for name, hits in hits.items():
        if len(hits) == 0: 
            continue
        ok = any(round(v[0][0], 2) >= 99.99 and round(v[0][3], 2) >= 99.99 for v in hits)
        if not ok: 
            partials[name] += hits

    tm = len(partials)
    print ':: Partial {} hits ({:.1f}%)'.format(tm, 100.0*tm/df.shape[0])
    for k in sorted(partials.keys(), key=lambda y: sum(yy[0][0] + yy[0][3] for yy in partials[y])):
        print '   -- partial http://humanparalogy.gs.washington.edu/build37/{0}'.format(k)
        for ((p1, n1, t1, p2, n2, t2), A, B) in sorted(partials[k]):
            print '      === {:.1f}% ({} of {}) and {:.1f}% ({} of {})'.format(p1, n1, t1, p2, n2, t2)
            A += (A[2] - A[1], A[5] - A[4])
            B += (B[2] - B[1], B[5] - B[4])
            print '          wgac:  {:5} {:11,}..{:11,} -> {:5} {:11,}..{:11,} ... len {:9,} -> {:9,} '.format(*A)
            print '          sedef: {:5} {:11,}..{:11,} -> {:5} {:11,}..{:11,} ... len {:9,} -> {:9,} '.format(*B)
except IOError:
    pass


# G:  37.0   2.8  40.8
# S:  19.7   2.7  22.4
# X:  74.0   2.9  77.1
# Y:  71.0   2.0  74.0
