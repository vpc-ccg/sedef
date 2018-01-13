/// 786

#include <bits/stdc++.h>
#include "common.h"
#include "fasta.h"
#include "align.h"
#include "search.h"

#include "extern/ksw2.h"

using namespace std;

// (hash, loci), sorted by loci
// vector<Minimizer> minimizers;
// hash --> list of locations 
// unordered_map<Hash, list<int>> index;

struct Anchor {
	list<int> query;
	list<int> ref;
	// list<int> hash;
};

// void chain_align(const string &sa, const string &sb, vector<Anchor> &chain)
// {
// 	int pa = chain[0].a, pb = chain[0].b;

// 	vector<pair<char, int>> cigar;
// 	for (int i = 1; i < chain.size(); i++) {
// 		Alignment tmp = align(sa.substr(pa, chain[0].a - pa), sb.substr(pb, chain[0].b - pb));
// 		cigar.insert(cigar.end(), tmp.cigar.begin(), tmp.cigar.end());
// 		cigar.push_back({'M', KMER_SIZE});
// 	}

// 	print >> cigar;
// }

auto chain(vector<Anchor> &v)
{
	auto similarity = [](int i, int j) {
		return 1;
	};

	vector<int> count(v.size(), 0);
	vector<int> prev(v.size(), 0);

	int argmax = 0;
	for (int i = 1; i < v.size(); i++) {
		for (int j = 0; j < i; j++) {
			int sim = similarity(i, j);
			if (count[j] + sim > count[i]) {
				count[i] = count[j] + sim;
				prev[i] = j;
			}
		}
		if (count[i] > count[argmax]) {
			argmax = i;
		}
	}

	vector<Anchor> result;
	do {
		result.push_back(v[argmax]);
		argmax = prev[argmax];
	} while (argmax != 0);

	return result;
}

auto cluster(const Index &hquery, const Index &href) 
{
	// MAP hb TO ha (hb=query; ha=ref)

	auto T = cur_time();

	vector<pair<int, int>> pairs; // sorted by query, then ref
	for (auto &query: hquery.minimizers) {
		auto ptr = href.index.find(query.hash);
		if (ptr == href.index.end()) 
			continue;

		for (auto ref_loc: ptr->second) { // TODO should be sorted --- check!
			pairs.push_back({query.loc, ref_loc});
		}
	}
	eprn("## pairs = {:n}", pairs.size());

	eprn(":: elapsed/pairs = {}s", elapsed(T)), T=cur_time();


	const int MAXGAP = 250;
	auto pless = [](const pair<int, int> &prev, const pair<int, int> &now) { // da li je [prev] < [now] ?
		if (prev.first + KMER_SIZE > now.first) 
			return false;
		if (prev.second + KMER_SIZE > now.second) 
			return false;
		if (now.first - prev.first > MAXGAP) 
			return false;
		if (now.second - prev.second > MAXGAP) 
			return false;
		return true;
	};

	for (int i = 1; i < pairs.size(); i++)
		assert(pairs[i - 1] < pairs[i]);

	int argmax = 0;
	vector<int> count(pairs.size(), 0), 
	            prev(pairs.size(), 0);
	
	for (int i = 1; i < pairs.size(); i++) {
		// for (int j = 0; j < i; j++) {
		for (int j = i - 1; j >= 0; j--) {
			if (pairs[i].first - pairs[j].first > MAXGAP) 
				break;
			bool ex = pless(pairs[j], pairs[i]);
			if (count[j] + 1 > count[i] && ex) {
				count[i] = count[j] + 1;
				prev[i] = j;
			}
		}
		if (count[i] > count[argmax]) {
			argmax = i;
		}
	}

	vector<list<int>> chains;
	vector<pair<int, int>> cntx(count.size());
	vector<char> used(count.size(), 0);
	for (int i = 0; i < count.size(); i++) 
		cntx[i] = {count[i], i};
	sort(cntx.begin(), cntx.end(), greater<pair<int,int>>());
	for (int i = 0; i < cntx.size(); i++) {
		if (cntx[i].first < 50)  // calculate this better!
			break;

		// Recover the path
		int start = cntx[i].second;
		if (used[start]) continue;
		list<int> chain;
		while (1) {
			// assert(!used[start]);
			used[start] = 1;
			chain.push_back(start);
			if (!start || !pless(pairs[prev[start]], pairs[start]) || used[prev[start]])
				break;
			start = prev[start];
		}
		if (chain.size() >= 50) {
			eprn("chain: len = {} / {}", chain.size(), cntx[i].first);
			eprn("   q: {}..{} <~> r: {}..{}", 
				pairs[chain.back()].first, pairs[chain.front()].first, 
				pairs[chain.back()].second, pairs[chain.front()].second
			);
			chains.push_back(chain);
		}
	}
	eprn(":: elapsed/chain = {}s [maxlen={}]", elapsed(T), count[argmax]), T=cur_time();

	for (auto &chain: chains) {
		int qs = pairs[chain.back()].first,  qe = pairs[chain.front()].first  + KMER_SIZE;
		int rs = pairs[chain.back()].second, re = pairs[chain.front()].second + KMER_SIZE;
		auto aln = align(hquery.seq->seq.substr(qs, qe - qs), href.seq->seq.substr(rs, re - rs),
			5,-4,40,1,max(qe-qs,re-rs)/4
		);
		auto err = aln.calculate_error();
		eprn("Aln: len={:n}, err={} (g={}, m={}), cigar={}", max(qe-qs,re-rs), 
			err.error(), err.gap_error(), err.mis_error(),
			aln.cigar_string());
	}
	eprn(":: elapsed/alignment = {}s", elapsed(T)), T=cur_time();

	// 3749 woo'hoo''

	// auto locate = [&](int i) {
	// 	// find all 
	// 	int lo = 0, hi = S.size() - 1;
	// 	while (lo <= hi) {
	// 		int mid = lo + (hi - lo) / 2;
	// 		if ()
	// 	}
	// };

	// vector<pair<int, int>> S { pairs[0] }; // Also sorted (query, ref)
	// count[0] = 1;
	// for (int i = 1; i < pairs.size(); i++) { 
	// 	if (pless(S.back(), pairs[i])) {
	// 		S.push_back(pairs[i]);
	// 	} else {
	// 		// returns first i s.t. X[i] >= y
	// 		auto pos = lower_bound(S.begin(), S.end(), pairs[i], pless);
	// 		*pos = pairs[i];
	// 	}

	// 	count[i] = S.size();
	// 	if (count[i] > count[argmax]) {
	// 		argmax = i;
	// 	}
	// }


	
	// clustering
	return true;
}


void test_align(const string &sa, const string &sb)
{
	Index ha(make_shared<Sequence>("A", sa)), 
	      hb(make_shared<Sequence>("B", sb));
	eprn("size = {} / {}", sa.size(), sb.size());

	auto T = cur_time();

	auto clusters = cluster(ha, hb);
	
	// for (auto &cluster: clusters) {
	// 	auto chain = chain(cluster);
	// 	auto aln = chain_align(chain);
	// }      

	eprn("elapsed={}s", elapsed(T));
}

void test(int, char** argv)
{
	FastaReference fr("data/hg19/chr2.fa");
	string s = "chr2 96307076 97264249 chr2 114040907 115040911";
	// ifstream fin(argv[0]);
	while (1) {
		auto ss = split(s, ' ');

		for (auto sss: ss) eprnn("{}_", sss); eprn("");

		const int OX = 500000;
		auto xa = fr.get_sequence(ss[0], atoi(ss[1].c_str())-OX, atoi(ss[2].c_str())+OX);
		auto xb = fr.get_sequence(ss[3], atoi(ss[4].c_str())-OX, atoi(ss[5].c_str())+OX);

		test_align(xa, xb);
		break;
	}
}

// ATAGCTAGCTAGCAT
// --AGCTAcC--GCATA