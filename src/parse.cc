/// 786

#include <bits/stdc++.h>
#include <seqan/seeds.h>

#include "common.h"
#include "fasta.h"
#include "align.h"
#include "search.h"

#include "extern/ksw2.h"

using namespace std;

/******************************************************************************/

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

// #define eprn(s,...) 0

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

/******************************************************************************/

auto find_chains(vector<pair<Anchor, bool>> &anchors)
{
	auto pless = [&](int i, int j) { // is anchors[i] < anchors[j] ?
		auto &pp = anchors[i].first;
		auto &p  = anchors[j].first;
		return pp.query + pp.query_len < p.query && pp.ref + pp.ref_len < p.ref;
	};

	vector<int> prev(anchors.size(), 0);
	vector<int> S { 0 };
	for (int i = 1; i < anchors.size(); i++) {
		if (pless(S.back(), i)) {
			prev[i] = S.back();
			S.push_back(i);
		} else {
			// returns first i s.t. X[i] >= y
			auto pos = distance(S.begin(), lower_bound(S.begin(), S.end(), i, pless));
			assert(pos > 0);
			prev[i] = S[pos - 1];
			S[pos] = i;
		}
	}

	eprn("-- len = {}", S.size());

	deque<int> path;
	vector<int> splits{0};

	int cur = S.back();
	for (int i = 0; i < S.size(); i++) {
		path.push_front(cur);
		cur = prev[cur];
	}


	auto pgap = [&](int i, int j) {
		auto &pp = anchors[i].first;
		auto &p  = anchors[j].first;

		const int MAXGAP = 250;

		int gap = p.query - (pp.query + pp.query_len);
		if (gap > MAXGAP && double(gap) / (p.query + p.query_len - pp.query) > MAX_ERROR) 
			return 2;
		gap = p.ref - (pp.ref + pp.ref_len);
		if (gap > MAXGAP && double(gap) / (p.ref + p.ref_len - pp.ref) > MAX_ERROR) 
			return 1;
		
		return 0;
	};
	for (int i = 1; i < path.size(); i++) {
		if (pgap(path[i - 1], path[i])) {
			splits.push_back(i);
		}
	}

	vector<vector<int>> paths;
	for (int i = 1; i < splits.size(); i++) {
		if (anchors[path[splits[i] - 1]].first.query - anchors[path[splits[i - 1]]].first.query + KMER_SIZE > 1000 && 
			anchors[path[splits[i] - 1]].first.ref   - anchors[path[splits[i - 1]]].first.ref   + KMER_SIZE > 1000) 
		{
			paths.push_back(vector<int>());
			for (int j = splits[i - 1]; j < splits[i]; j++) {
				paths.back().push_back(path[j]);
				anchors[path[j]].second = false;
			}
		}
	}

	eprn("-- chains = {}", paths.size());
	return paths;
}

auto cluster(const Index &hquery, const Index &href) 
{
	auto T = cur_time();

	vector<pair<Anchor, bool>> pairs; // sorted by query, then ref

	eprn("-- clustering len {:n} vs {:n}", hquery.seq->seq.size(), href.seq->seq.size());
	for (auto &query: hquery.minimizers) {
		auto ptr = href.index.find(query.hash);
		if (ptr == href.index.end()) 
			continue;
		for (auto ref_loc: ptr->second) {
			pairs.push_back({{
				query.loc, ref_loc, 
				hquery.kmer_size, href.kmer_size, 
				{query.loc}, {ref_loc}
			}, true});
			// addSeed(seeds, Seed<Simple>(query.loc, query.loc+KMER_SIZE, ref_loc, ref_loc+KMER_SIZE), Single());
			if (pairs.size() > 1) {
				assert(pairs[pairs.size() - 2].first < pairs.back().first);
			}
		}
	}

	vector<Anchor> paths;
	for (int iter = 0; iter < 10; iter++) {
		eprn("-- pairs = {:n}", pairs.size());
		auto pp = find_chains(pairs);
		if (pp.size() == 0) 
			break;

		for (auto &path: pp) {
			Anchor a;

			int qlo = pairs[path.front()].first.query, qhi = pairs[path.back()].first.query;
			int rlo = pairs[path.front()].first.ref,   rhi = pairs[path.back()].first.ref;
			// eprn("{:n}..{:n} --- {:n}..{:n}", qlo, qhi, rlo, rhi);

			a.query = qlo;
			a.query_len = qhi - qlo + KMER_SIZE;
			a.ref = rlo;
			a.ref_len = rhi - rlo + KMER_SIZE;
			for (int pi: path) {
				a.query_kmers.push_back(pairs[pi].first.query_kmers.front());
				a.ref_kmers.push_back(pairs[pi].first.ref_kmers.front());
			}
			for (int pi = path.front(); pi <= path.back(); pi++) {
				assert(pairs[pi].first.query >= qlo);
				assert(pairs[pi].first.query <= qhi);
				if (pairs[pi].first.ref >= rlo && pairs[pi].first.ref <= rhi) {
					pairs[pi].second = false;
				}
			}
			paths.push_back(a);
		}

		vector<pair<Anchor, bool>> px;
		for (auto &p: pairs) 
			if (p.second) px.push_back(p);
		pairs = px;
		eprn(":: iter = {}, elapsed/dp = {}s", iter, elapsed(T)), T=cur_time();
	}

	// eprn(":: elapsed/pairs = {}s", elapsed(T)), T=cur_time();
	// auto pless = [&](const Anchor &pp, const Anchor &p) {
	// 	// Ignore overlapping pairs
	// 	if (pp.query + pp.query_len > p.query) 
	// 		return 1;
	// 	if (pp.ref + pp.ref_len > p.ref) 
	// 		return 1;

	// 	// Make sure to respect gap penalty
	// 	const int MAXGAP = 250;
	// 	int gap = p.query - (pp.query + pp.query_len);
	// 	if (gap > MAXGAP && double(gap) / (p.query + p.query_len - pp.query) > MAX_ERROR) 
	// 		return 2;
 
	// 	gap = p.ref - (pp.ref + pp.ref_len);
	// 	if (gap > MAXGAP && double(gap) / (p.ref + p.ref_len - pp.ref) > MAX_ERROR) 
	// 		return 1;
		
	// 	return 0;
	// };

	// // Run DP
	// vector<pair<int, int>> count(pairs.size(), {0, 0}); 
	// vector<int> prev(pairs.size(), 0);
	// for (int i = 1; i < pairs.size(); i++) {
	// 	count[i].second = i;
	// 	for (int j = i - 1; j >= 0; j--) {
	// 		char status = pless(pairs[j], pairs[i]);
	// 		if (status == 1) 
	// 			continue;
	// 		else if (status == 2) 
	// 			break;
	// 		// Update DP table
	// 		if (count[j].first + 1 > count[i].first) {
	// 			count[i].first = count[j].first + 1;
	// 			prev[i] = j;
	// 		}
	// 	}
	// }
	// eprn(":: elapsed/dp = {}s", elapsed(T)), T=cur_time();

	// // Reconstruct long chains
	// vector<Anchor> new_pairs;
	// vector<bool> used(count.size(), 0);
	// sort(count.begin(), count.end(), greater<pair<int,int>>());
	// for (auto &c: count) { // From highest to the lowest
	// 	// if (c.first < 50)  // TODO calculate this better!
	// 		// break;
	// 	// Recover the path
	// 	int start = c.second;
	// 	if (used[start]) 
	// 		continue;
	// 	Anchor extension = pairs[start];
	// 	used[start] = true;
	// 	while (start && !used[prev[start]] && !pless(pairs[prev[start]], pairs[start])) {
	// 		start = prev[start];
	// 		used[start] = true;

	// 		extension.query_len += extension.query - pairs[start].query;
	// 		extension.query = pairs[start].query;

	// 		extension.ref_len += extension.ref - pairs[start].ref;
	// 		extension.ref = pairs[start].ref;

	// 		extension.query_kmers.insert(extension.query_kmers.begin(), pairs[start].query_kmers.begin(), pairs[start].query_kmers.end());
	// 		extension.ref_kmers.insert(extension.ref_kmers.begin(), pairs[start].ref_kmers.begin(), pairs[start].ref_kmers.end());
	// 	}
	// 	if (extension.query_kmers.size() > 1 
	// 		&& extension.query_len >= 1000 
	// 		&& extension.ref_len >= 1000) 
	// 	{
	// 		new_pairs.push_back(extension);
	// 	}
	// }
	// eprn("-- found {} chains", new_pairs.size());
	// eprn(":: elapsed/chain = {}s", elapsed(T)), T=cur_time();


	vector<Hit> hits;
	for (auto &p: paths) {
		auto aln = p.anchor_align(hquery.seq->seq, href.seq->seq);
		hits.push_back({
			hquery.seq, p.query, p.query + p.query_len,
			href.seq, p.ref, p.ref + p.ref_len,
			9999, "", "",
			aln
		});
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
	}
	eprn(":: elapsed/alignment = {}s", elapsed(T)), T=cur_time();
	
	return hits;
}

/******************************************************************************/

vector<Hit> test_align(const string &sa, const string &sb)
{
	return cluster(
		make_shared<Sequence>("A", sa), 
		make_shared<Sequence>("B", sb)
	);
}

void test(int, char** argv)
{
	FastaReference fr("data/hg19/hg19.fa");
	// 200sec
	string s = "chr22	20318686	20514796	chr22	21513387	21793546		-1	+	+	280159	0		2025	OK";
	// ifstream fin(argv[0]);
	while (1) {
		Hit h = Hit::from_bed(s);
		
		auto T = cur_time();
		const int OX = 0;
		auto q = fr.get_sequence(h.query->name, h.query_start, h.query_end);
		auto r = fr.get_sequence(h.ref->name, h.ref_start, h.ref_end);
		test_align(q, r);
		eprn("total {}s", elapsed(T));
		break;
	}
}

// ATAGCTAGCTAGCAT
// --AGCTAcC--GCATA