# 786

#%%
import pysam, os, yaml, re, copy, subprocess, tempfile, glob, sys
import seaborn as sns
from collections import *
import cPickle as pickle
import pandas as pd
import numpy as np
import subprocess
import scipy
from collections import defaultdict
import rpy2
%load_ext rpy2.ipython
%matplotlib inline
def maybeint(s):
    try: return int(s)
    except ValueError: return s
pd.set_option("display.max_rows",10)
pd.set_option("display.max_columns",100)


#%% Calculate the scores
%%bash
g++ -o newidea newidea.cc -O2 -std=c++14 patterns.o fmt/fmt/format.o
cat data/GRCh37GenomicSuperDup.tab | awk '{print $1,$2,$3,$17,$5,"+",$7,$8,$9,$4,$5,$6}' 'OFS=\t' \
    | grep -v '_gl' | tr '_' '-' | sed 1d | sed 's/\([+-]\)\t/\1\n/' \
    | bedtools getfasta -name -fi data/hg19/hg19.fa -bed - -fo - -s  | head -n4 | ./newidea > newidea_data.txt

#%% Load data
df = pd.read_table("sedef/data/GRCh37GenomicSuperDup.tab")
dscores = pd.read_table("sedef/newidea_data.txt", names=['alignfile', 'name', 'coreSet', 'qgramSet', 'spaced1Set', 'spaced2Set'], index_col=False)
dscores['alignfile'] = dscores.alignfile.apply(lambda x: x[1:].replace('-', '_'))
dscores['name'] = dscores.name.apply(lambda x: x[1:])
df = pd.merge(df, dscores, on=['alignfile', 'name'])
df.head(2)

#%% Connect WGAC data with alignments *AND* masked regions
def get_seq(r1, r2, s1, s2):
    r1, r2 = map(lambda x: x.replace('\t', '\\t'), [r1, r2])
    q = subprocess.check_output('echo "{}\n{}" | bedtools getfasta -fi sedef/data/hg19/hg19.fa -bed - -fo - -s'.format(r1, r2), shell=True)
    q = map(str.strip, q.split())
    def apply_lower(S):
        s, r = S
        ib = 0
        s = list(s)
        for ic in xrange(len(s)):
            if s[ic] == '-': continue
            s[ic] = r[ib]
            ib += 1
        return ''.join(s)
    return map(apply_lower, [(s1, q[1]), (s2, q[3])])
seqs = [[], []]
for ri, r in df.iterrows():
    af = r.alignfile.split('/')[2]
    with open("sedef/alignments/" + af) as f:
        fl = [l.strip().split()[1:] for l in f.readlines() if l[0] in 'cF']

    rs, re = int(fl[0][2][1:]), int(fl[0][4][:-1])
    r1 = '{}\t{}\t{}\tA\t0\t+'.format(fl[0][1], rs - 1, re)
    rs, re, st2 = int(fl[0][7][1:]), int(fl[0][9][:-1]), '+'
    if rs > re: rs, re, st2 = re + 1, rs + 1, '-'
    r2 = '{}\t{}\t{}\tB\t0\t{}'.format(fl[0][6], rs - 1, re, st2)
    seq1 = ''.join(fl[i][0] for i in range(1, len(fl), 2))
    seq2 = ''.join(fl[i][0] for i in range(2, len(fl), 2))
    seq1, seq2 = get_seq(r1, r2, seq1, seq2)
    seqs[0].append(seq1)
    seqs[1].append(seq2)
df = df.assign(seq1 = seqs[0]).assign(seq2 = seqs[1])
df.to_csv("sedef/data/GRCh37GenomicSuperDup_Extra.csv", sep="\t")


#%%
data = [[], []]
for ri, r in df.iterrows():
    for i in xrange():
    norepe[0].append(sum(1 for c in r.seq1 if c.islower()))
    norepe[1].append(sum(1 for c in r.seq1 if c.islower()))
df = df.assign(seq1 = seqs[0]).assign(seq2 = seqs[1])

#%% Plot the data
# Plot the cores / 1000bp!
sns.distplot(df.coreSet / (df.otherSize / 750), kde=False)


dp = df[df.coreSet / (df.otherSize / 750) < 15]
sns.distplot(dp.matchB, kde=False)
100.0 * dp.shape[0] / df.shape[0] # < 3%

sns.distplot(df.qgramSet / (df.otherSize / 1000), kde=False)
sns.distplot(df.spaced2Set / (df.otherSize / 1000), kde=False)

df

#%% Parse the data
stats = {}

for fn in sorted(glob.glob('../sedef/data/wgac/build37/align_both/0001/both*')):
    name = os.path.basename(fn)
    print name
    sa, sb = [], []
    soa, sob = '', ''
    with open(fn) as f:
        l = next(f).split()[2:]
        cA, sA, eA = l[0], l[1][1:], l[3][:-1]
        l = l[5:]
        cB, sB, eB = l[0], l[1][1:], l[3][:-1]
        next(f)
        for li, l in enumerate(f):
            if li % 6 == 2: sa.append(l.split()[1])
            if li % 6 == 4: sb.append(l.split()[1])
    if 'random' in cA or 'random' in cB: continue
    if 'chrUn' in cA or 'chrUn' in cB: continue
    soa = get_seq(cA, sA, eA)
    sob = get_seq(cB, sB, eB)
    sa = list(''.join(sa))
    sb = list(''.join(sb))
    assert(len(sa) == len(sb))
    gap_len = 0
    mut_len = 0
    free_len = 0
    gaps = []
    frees = []
    for ca, cb in zip(sa, sb):
        if ca != cb:
            if free_len > 0: frees.append(free_len)
            free_len = 0
        if ca == '-' or cb == '-':
            gap_len += 1
        else:
            if gap_len != 0: gaps.append(gap_len)
            gap_len = 0
            if ca != cb:
                mut_len += 1
            else:
                free_len += 1
    ia, ib = 0, 0 # indices in "real" seqs
    gap_len_fug = 0
    mut_len_fug = 0
    free_len_fug = 0
    gaps_fug = []
    frees_fug = []
    # print name
    assert(''.join(c for c in sa if c != '-') == soa.upper())
    assert(''.join(c for c in sb if c != '-') == sob.upper())
    i = 0
    for ci, c in enumerate(sa):
        if c != '-' and soa[i].islower(): sa[ci] = '-' # c.lower()
        if c != '-': i += 1
    i = 0
    for ci, c in enumerate(sb):
        if c != '-' and sob[i].islower(): sb[ci] = '-' # c.lower()
        if c != '-': i += 1
    nl = 0
    for ca, cb in zip(sa, sb):
        if ca == '-' and cb == '-':
            continue
        nl += 1
        if ca != cb:
            if free_len_fug > 0:
                frees_fug.append(free_len_fug)
            free_len_fug = 0
        if ca == '-' or cb == '-':
            gap_len_fug += 1
        else:
            if gap_len_fug != 0: gaps_fug.append(gap_len_fug)
            gap_len_fug = 0
            if ca != cb:
                mut_len_fug += 1
            else:
                free_len_fug += 1

    stats[name] = (len(sa), mut_len, gaps, frees, nl, mut_len_fug, gaps_fug, frees_fug)

#%% just load it... easier
pickle.dump(stats, open('wgac-chr1-stats.pickle', 'w'))

#%% Calculate stats
mut_rates = []   # how many mutations per 1000 "normal" bps
gap_perc = []    # how many gap bases per 1000 bps
gap_counts = []  # how many gaps per 1000 "normal" bps
free_lens = []   # length of non-error block
gap_lens = []    # length of gap block
f_mut_rates = []   # how many mutations per 1000 "normal" bps
f_gap_perc = []    # how many gap bases per 1000 bps
f_gap_counts = []  # how many gaps per 1000 "normal" bps
f_free_lens = []   # length of non-error block
f_gap_lens = []    # length of gap block
f_lens = []
for _, (l, ml, gaps, frees, f_l, f_ml, f_gaps, f_frees) in stats.iteritems():
    #print _, ml, l - sum(gaps)
    mut_rates.append(float(ml) / (float(l - sum(gaps)) / 1000))
    gap_counts.append(float(len(gaps)) / (float(l - sum(gaps)) / 1000))
    gap_perc.append(float(sum(gaps)) / (float(l) / 1000))
    gap_lens += gaps
    free_lens += [frees]
    if f_l > 0:
        print _, f_ml, f_l, sum(f_gaps)
        f_mut_rates.append(float(f_ml) / (float(f_l - sum(f_gaps)) / 1000))
        f_gap_counts.append(float(len(f_gaps)) / (float(f_l - sum(f_gaps)) / 1000))
        f_gap_perc.append(float(sum(f_gaps)) / (float(f_l) / 1000))
        f_free_lens += f_frees
        f_gap_lens += f_gaps
        f_lens.append(f_l)
mut_rates = np.array(mut_rates)
gap_perc = np.array(gap_perc)
gap_counts = np.array(gap_counts)
free_lens = np.array(free_lens)
gap_lens = np.array(gap_lens)
f_mut_rates = np.array(f_mut_rates)
f_gap_perc = np.array(f_gap_perc)
f_gap_counts = np.array(f_gap_counts)
f_free_lens = np.array(f_free_lens)
f_gap_lens = np.array(f_gap_lens)
f_lens = np.array(f_lens)

#%% Plot percentage of mutations
sns.distplot(100 * (mut_rates / 1000.0))
sns.distplot(100 * (f_mut_rates / 1000.0))
# dafuq case: both005477


#%%
g = sns.distplot(gap_counts, bins=1000)
g.set(xlim=(-1, 25))
g
sns.distplot(f_gap_counts, bins=1000)


#%%
g = sns.distplot([w for w in gap_perc if w > 0], bins=1000)
g.set(xlim=(-1, 200))
g
sns.distplot([w for w in f_gap_perc if w > 0], bins=1000)

#%%
g = sns.distplot(free_lens, bins=2000)
g.set(xlim=(-1, 50))
g
sns.distplot(f_free_lens, bins=2000)

#%%
Counter(free_lens)
