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
        A = [r.chrom, r.chromStart, r.chromEnd,
            r.otherChrom, r.otherStart, r.otherEnd,
            r.alignfile, 0, '+', '-' if r.strand == '_' else '+']
        print >>f, '\t'.join(map(str, A))
        hits[r.alignfile] = list()
        name_to_coor[r.alignfile] = r #(r.chromStart, r.chromEnd, r.otherStart, r.otherEnd)


def overlap(sa, ea, sb, eb):
    return max(0, min(ea, eb) - max(sa, sb))

def match(sW, eW, sS, eS):
    oo = overlap(sW, eW, sS, eS)
    wW = oo
    dW = 100.0 * wW / float(eW - sW)
    nE = (eS-sS+eW-sW-wW) / float(eS-sS) # how much extention is needed in pct of sedef's size
    return [(dW, wW, eW - sW, nE)]

# print overlap(60426,69533, 61131,    111437)
# exit(0)

def diff(wgac, sedef): # how much wgac is off sedef
    if not tuple(wgac[0:3]) < tuple(wgac[3:6]):
        wgac = wgac[3:6] + wgac[0:3] + wgac[6:]
    if not tuple(sedef[0:3]) < tuple(sedef[3:6]):
        sedef = sedef[3:6] + sedef[0:3] + sedef[6:]
    return match(*(wgac[1:3] + sedef[1:3])) + match(*(wgac[4:6] + sedef[4:6]))

def process_path(path, pnew):
    with open(pnew, 'w') as fw:
        with open(path) as f:
            for l in f:
                l = l.strip().split('\t')
                s1, e1, s2, e2 = map(int, l[1:3] + l[4:6])
                o = max(e1-s1,e2-s2)*1
                l[1:3] = [max(1, s1 - o), e1 + o]
                l[4:6] = [max(1, s2 - o), e2 + o]
                print >>fw, '\t'.join(map(str, l))

pnew = path #+ "____"
# process_path(path, pnew)
system("bedtools pairtopair -a temp.bed -b <(cat {} | tr -d ,) -type both > temp_diff.bed".format(pnew))
# os.unlink(pnew)
print ':: After bedtools we have {} hits to process'.format(system("wc -l temp_diff.bed"))

with open('temp_diff.bed') as f:
    for l in f:
        l = l.strip().split()
        q = 0
        A = (l[q], int(l[q+1]), int(l[q+2]), l[q+3], int(l[q+4]), int(l[q+5]), l[q+8]+l[q+9]) # WGAC
        q = 10
        # print [(i,x) for i,x in enumerate(l)]
        B = (l[q], int(l[q+1]), int(l[q+2]), l[q+3], int(l[q+4]), int(l[q+5]), l[q+6]+l[q+7]) # SEDEF

        # B = list(B)
        # d = min(20000, 4 * max(B[2] - B[1], B[5] - B[4]))
        # B[1] -= d
        # B[2] += d
        # B[4] -= d
        # B[5] += d
        # B = tuple(B)
        # exit(0)

        name = l[6]
        hits[name].append((diff(A, B), A, B))

try:
    tm = sum(1 for k, vs in hits.items() if len(vs) == 0)
    print ':: Missed {} hits ({:.1f}%)'.format(tm, 100.0*tm/len(hits))
    for name, vs in sorted(hits.items()):
        if len(vs) == 0:
            print '   -- missed http://humanparalogy.gs.washington.edu/build37/{0}'.format(name)
            r = name_to_coor[name]
            print '      wgac: {:5} {:11,}..{:11,} -> {:5} {:11,}..{:11,} ... {} ... len {:9,} -> {:9,} '.format(
                r.chrom, r.chromStart, r.chromEnd, r.otherChrom, r.otherStart, r.otherEnd, r.strand,
                r.chromEnd - r.chromStart, r.otherEnd - r.otherStart
            )

    partials = defaultdict(list)
    for name, h in hits.items():
        if len(h) == 0: 
            continue
        ok = any(v[0][0][0] >= 100 and v[0][1][0] >= 100 for v in h)
        if not ok: 
            partials[name] += h

    tm = len(partials)
    print ':: Partial {} hits ({:.1f}%)'.format(tm, 100.0*tm/len(hits))
    for k in sorted(partials.keys(), key=lambda y: sum(yy[0][0][0] + yy[0][1][0] for yy in partials[y])):
        print '   -- partial http://humanparalogy.gs.washington.edu/build37/{0}'.format(k)
        for (((p1, n1, t1, e1), (p2, n2, t2, e2)), A, B) in sorted(partials[k]):
            print '      === {:.1f}% ({} of {}) and {:.1f}% ({} of {})'.format(p1, n1, t1, p2, n2, t2)
            A += (A[2] - A[1], A[5] - A[4])
            B += (B[2] - B[1], B[5] - B[4])
            print '          wgac:  {:5} {:11,}..{:11,} -> {:5} {:11,}..{:11,} {} ... len {:9,} -> {:9,} '.format(*A)
            print '          sedef: {:5} {:11,}..{:11,} -> {:5} {:11,}..{:11,} {} ... len {:9,} -> {:9,} '.format(*B)
            # print '>> need: {:.2f} {:.2f}'.format(e1, e2)

    tm = sum(1 for k, v in hits.iteritems() if k not in partials and len(v) > 0)
    print ':: Full {} hits ({:.1f}%)'.format(tm, 100.0*tm/len(hits))
    for k, v in hits.iteritems():
        if k in partials or len(v) == 0:
            continue

        found = 0
        for (((p1, n1, t1, e1), (p2, n2, t2, e2)), A, B) in sorted(v, reverse=True):
            if not tuple(A[0:3]) < tuple(A[3:6]):
                A = A[3:6] + A[0:3] + (A[6],)
            if not tuple(B[0:3]) < tuple(B[3:6]):
                B = B[3:6] + B[0:3] + (B[6],)

            chr_eq = (A[0] == B[0] and A[3] == B[3])
            str_eq = A[6] == B[6]
            l_eq = A[1] >= B[1] and A[2] <= B[2] 
            r_eq = A[4] >= B[4] and A[5] <= B[5]

            eq = chr_eq and str_eq and l_eq and r_eq
            if eq:
                found = 1
                break
        if found == 0:
            print 'whooops!!!'
            for (((p1, n1, t1, e1), (p2, n2, t2, e2)), A, B) in sorted(v, reverse=True):
                print p1, p2
                print n1, n2
                print t1, t2
                print '          wgac:  {:5} {:11,}..{:11,} -> {:5} {:11,}..{:11,} {} '.format(*A)
                print '          sedef: {:5} {:11,}..{:11,} -> {:5} {:11,}..{:11,} {} '.format(*B)
            exit(0)


except IOError:
    pass


# G:  37.0   2.8  40.8
# S:  19.7   2.7  22.4
# X:  74.0   2.9  77.1
# Y:  71.0   2.0  74.0
