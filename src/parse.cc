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
	int query;
	int ref;
	int query_len;
	int ref_len;
	list<int> query_kmers;
	list<int> ref_kmers;

	bool operator< (const Anchor &a) const {
		return tie(query, ref) < tie(a.query, a.ref);
 	}
 	Alignment anchor_align(const string &qs, const string &rs) ;
};

void append_cigar(Alignment &src, const deque<pair<char, int>> &app)
{
	assert(app.size());
	assert(src.cigar.size());
	if (src.cigar.back().first == app.front().first) {
		src.cigar.back().second += app.front().second;
		src.cigar.insert(src.cigar.end(), next(app.begin()), app.end());
	} else {
		src.cigar.insert(src.cigar.end(), app.begin(), app.end());
	}
}

Alignment Anchor::anchor_align(const string &qs, const string &rs) 
{
	const int kmer_size = KMER_SIZE;

	auto aln_from_match = [&](int qk, int rk) {
		return Alignment {
			"A", qk, qk + kmer_size, 
			"B", rk, rk + kmer_size,
			qs.substr(qk, kmer_size),
			rs.substr(rk, kmer_size),
			"", "", "",
			{{'M', kmer_size}}, {}
		};
	};

	auto qk_p = query_kmers.begin();
	auto rk_p = ref_kmers.begin();
	Alignment aln = aln_from_match(*qk_p, *rk_p);
	assert(query_kmers.size() == ref_kmers.size());
	for (auto qk = next(qk_p), rk = next(rk_p); qk != query_kmers.end(); qk++, rk++) {
		aln.end_a = *qk + kmer_size;
		aln.end_b = *rk + kmer_size;
		aln.a += qs.substr(*qk_p + kmer_size, *qk - *qk_p);
		aln.b += rs.substr(*rk_p + kmer_size, *rk - *rk_p);
		if (*qk - *qk_p - kmer_size && *rk - *rk_p - kmer_size) {
			auto gap = align(
				qs.substr(*qk_p + kmer_size, *qk - *qk_p - kmer_size), 
				rs.substr(*rk_p + kmer_size, *rk - *rk_p - kmer_size), 
				5, -4, 40, 1
			);
			append_cigar(aln, gap.cigar);
		} else if (*qk - *qk_p - kmer_size) {
			append_cigar(aln, {{'D', *qk - *qk_p - kmer_size}});	
		} else if (*rk - *rk_p - kmer_size) {
			append_cigar(aln, {{'I', *rk - *rk_p - kmer_size}});	
		} 
		append_cigar(aln, {{'M', kmer_size}});
		qk_p = qk, rk_p = rk;
	}

	aln.populate_nice_alignment();
	return aln;
}

auto cluster(const Index &hquery, const Index &href) 
{
	// MAP hb TO ha (hb=query; ha=ref)

	auto T = cur_time();

	vector<Anchor> pairs; // sorted by query, then ref
	for (auto &query: hquery.minimizers) {
		auto ptr = href.index.find(query.hash);
		if (ptr == href.index.end()) 
			continue;
		for (auto ref_loc: ptr->second) {
			pairs.push_back({query.loc, ref_loc, 
				hquery.kmer_size, href.kmer_size, 
				{query.loc}, {ref_loc}});
			if (pairs.size() > 1) {
				assert(pairs[pairs.size() - 2] < pairs.back());
			}
		}
	}
	eprn("-- pairs = {:n}", pairs.size());

	eprn(":: elapsed/pairs = {}s", elapsed(T)), T=cur_time();

	auto pless = [&](const Anchor &pp, const Anchor &p) {
		// Ignore overlapping pairs
		if (pp.query + pp.query_len > p.query) 
			return 1;
		if (pp.ref + pp.ref_len > p.ref) 
			return 1;

		// Make sure to respect gap penalty
		const int MAXGAP = 250;
		int gap = p.query - (pp.query + pp.query_len);
		if (gap > MAXGAP && double(gap) / (p.query + p.query_len - pp.query) > MAX_ERROR) 
			return 2;
 
		gap = p.ref - (pp.ref + pp.ref_len);
		if (gap > MAXGAP && double(gap) / (p.ref + p.ref_len - pp.ref) > MAX_ERROR) 
			return 1;
		
		return 0;
	};

	auto Ti = cur_time();
	for (int iter = 0; iter < 2; iter++) {	
	// Run DP
		vector<pair<int, int>> count(pairs.size(), {0, 0}); 
		vector<int> prev(pairs.size(), 0);
		for (int i = 1; i < pairs.size(); i++) {
			count[i].second = i;
			for (int j = i - 1; j >= 0; j--) {
				char status = pless(pairs[j], pairs[i]);
				if (status == 1) 
					continue;
				else if (status == 2) 
					break;
				// Update DP table
				if (count[j].first + 1 > count[i].first) {
					count[i].first = count[j].first + 1;
					prev[i] = j;
				}
			}
		}
		eprn("   :: elapsed/dp = {}s", elapsed(T)), T=cur_time();

	// Reconstruct long chains
		vector<Anchor> new_pairs;
		vector<bool> used(count.size(), 0);
		sort(count.begin(), count.end(), greater<pair<int,int>>());
		for (auto &c: count) { // From highest to the lowest
			// if (c.first < 50)  // TODO calculate this better!
				// break;

			// Recover the path
			int start = c.second;
			if (used[start]) 
				continue;
			Anchor extension = pairs[start];
			used[start] = true;
			while (start && !used[prev[start]] && !pless(pairs[prev[start]], pairs[start])) {
				start = prev[start];
				used[start] = true;

				extension.query_len += extension.query - pairs[start].query;
				extension.query = pairs[start].query;

				extension.ref_len += extension.ref - pairs[start].ref;
				extension.ref = pairs[start].ref;

				extension.query_kmers.insert(extension.query_kmers.begin(), pairs[start].query_kmers.begin(), pairs[start].query_kmers.end());
				extension.ref_kmers.insert(extension.ref_kmers.begin(), pairs[start].ref_kmers.begin(), pairs[start].ref_kmers.end());
			}
			if (extension.query_kmers.size() > 1)
				new_pairs.push_back(extension);
		}
		eprn("   -- found {} chains", new_pairs.size());
		eprn("   :: elapsed/chain = {}s", elapsed(T)), T=cur_time();

		pairs = new_pairs;
		eprn(":: elapsed/iter{:02d} = {}s", iter, elapsed(Ti)), Ti=cur_time();
	}

	eprn("-- total {} chains", pairs.size());
	// for (auto &p: pairs) {
	// 	if (p.query_kmers.size() < 100) continue;
	// 	auto aln = align(
	// 		hquery.seq->seq.substr(p.query, p.query_len), 
	// 		href.seq->seq.substr(p.ref, p.ref_len),
	// 		5, -4, 40, 1, 
	// 		max(p.query_len, p.ref_len) / 4
	// 	);
	// 	auto err = aln.calculate_error();
	// 	eprn("      ||> L {:n} / {:n} \n"
	// 		 "          Q {:n}..{:n} <~> R {:n}..{:n} \n"
	// 		 "          E {:4.2f} (g={:4.2f}, m={:4.2f})", 
	// 		p.query_len, p.ref_len,
	// 		p.query, p.query + p.query_len,
	// 		p.ref, p.ref + p.ref_len,
	// 		err.error(), err.gap_error(), err.mis_error()
	// 	);
	// }
	// eprn(":: elapsed/alignment = {}s", elapsed(T)), T=cur_time();

	for (auto &p: pairs) {
		if (p.query_kmers.size() < 100) continue;
		auto aln = p.anchor_align(hquery.seq->seq, href.seq->seq);
		auto err = aln.calculate_error();
		eprn("      ||> L {:n} / {:n} \n"
			 "          Q {:n}..{:n} <~> R {:n}..{:n} \n"
			 "          E {:4.2f} (g={:4.2f}, m={:4.2f}) --- {}", 
			p.query_len, p.ref_len,
			p.query, p.query + p.query_len,
			p.ref, p.ref + p.ref_len,
			err.error(), err.gap_error(), err.mis_error(),
			p.query_kmers.size() * KMER_SIZE
		);
		// eprn("{}", aln.print(60));
	}
	eprn(":: elapsed/alignment = {}s", elapsed(T)), T=cur_time();
	
	return true;
}


void test_align(const string &sa, const string &sb)
{
	Index ha(make_shared<Sequence>("A", sa)), 
	      hb(make_shared<Sequence>("B", sb));
	eprn("|| size = {:n} ~ {:n}", sa.size(), sb.size());

	auto T = cur_time();

	auto clusters = cluster(ha, hb);
	
	// for (auto &cluster: clusters) {
	// 	auto chain = chain(cluster);
	// 	auto aln = chain_align(chain);
	// }      

	eprn("|> elapsed/total = {}s", elapsed(T));
}

void test(int, char** argv)
{
	FastaReference fr("data/hg19/chr2.fa");
	string s = "chr2 96307076 97264249 chr2 114040907 115040911";
	// ifstream fin(argv[0]);
	while (1) {
		auto ss = split(s, ' ');

		for (auto sss: ss) eprnn("{} ", sss);

		const int OX = 0;
		auto xa = fr.get_sequence(ss[0], atoi(ss[1].c_str())-OX, atoi(ss[2].c_str())+OX);
		auto xb = fr.get_sequence(ss[3], atoi(ss[4].c_str())-OX, atoi(ss[5].c_str())+OX);

		test_align(xa, xb);
		break;
	}
}

// ATAGCTAGCTAGCAT
// --AGCTAcC--GCATA