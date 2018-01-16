/// 786

#include <bits/stdc++.h>
#include <seqan/seeds.h>

#include "common.h"
#include "fasta.h"
#include "align.h"
#include "search.h"

using namespace std;

/******************************************************************************/

struct Anchor {
	int query_start, query_end;
	int ref_start, ref_end;
	list<int> query_kmers, ref_kmers;

	bool operator< (const Anchor &a) const {
		return tie(query_start, ref_start) < tie(a.query_start, a.ref_start);
	}
	Alignment anchor_align(const string &qs, const string &rs) ;
};

#define eprn(s,...) 0

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
	aln.error = aln.calculate_error();
	return aln;
}

/******************************************************************************/

auto merge_long_chains(vector<Anchor> &anchors, const Index &query_hash, const Index &ref_hash) 
{
	assert(anchors.size());
	sort(anchors.begin(), anchors.end());

	vector<Alignment> alignments;
	for (auto &anchor: anchors) {
		alignments.push_back(anchor.anchor_align(query_hash.seq->seq, ref_hash.seq->seq));
	}

	vector<int> score, dp; // fill the scores!!
	for (int i = 0; i < anchors.size(); i++) {
		auto &err = alignments[i].error;
		score.push_back(err.matches * 5 - err.mismatches * 4 - 40 * err.gaps - err.gap_bases);
		dp.push_back(score[i]);
	}
	vector<int> prev(anchors.size(), -1);
	// maximize span while keeping the properties
	for (int i = 0; i < anchors.size(); i++) {
		for (int j = 0; j < i; j++) {
			if (anchors[j].query_end >= anchors[i].query_start) 
				continue;
			if (anchors[j].ref_end >= anchors[i].ref_start) 
				continue;

			int h = anchors[i].query_start - anchors[j].query_end;
			int v = anchors[i].ref_start - anchors[j].ref_end;
			int dist = score[i] - (4 * min(h, v) + 40 + abs(h - v));

			if (dp[j] + dist > dp[i]) {
				dp[i] = dp[j] + dist;
				prev[i] = j;
			}
		}
	}

	vector<int> dpidx(dp.size());
	for (int i = 0; i < dp.size(); i++) 
		dpidx[i] = i;
	sort(dpidx.begin(), dpidx.end(), [&](int i, int j) { return dp[i] > dp[j]; });
	
	vector<Hit> hits;
	for (int cur: dpidx) {
		if (dp[cur] < 0) 
			continue;

		// reconstruct
		deque<int> dq;
		while (cur != -1 && dp[cur] != -1) {
			dq.push_front(cur);
			dp[cur] = -1;
			cur = prev[cur];
		}
		assert(dq.size());
		
		Alignment aln = alignments[dq[0]];
		auto prev = dq.begin();
		for (auto cur = next(prev); cur != dq.end(); cur++) {
			aln.end_a = anchors[*cur].query_end;
			aln.end_b = anchors[*cur].ref_end;
			aln.a += query_hash.seq->seq.substr(anchors[*prev].query_end, anchors[*cur].query_end - anchors[*prev].query_end);
			aln.b += ref_hash.seq->seq.substr(anchors[*prev].ref_end, anchors[*cur].ref_end - anchors[*prev].ref_end);

			if (anchors[*cur].query_start - anchors[*prev].query_end && anchors[*cur].ref_start - anchors[*prev].ref_end) {
				auto gap = align(
					query_hash.seq->seq.substr(anchors[*prev].query_end, anchors[*cur].query_start - anchors[*prev].query_end),
					ref_hash.seq->seq.substr(anchors[*prev].ref_end, anchors[*cur].ref_start - anchors[*prev].ref_end),
					5, -4, 40, 1
				);
				append_cigar(aln, gap.cigar);
			} else if (anchors[*cur].query_start - anchors[*prev].query_end) {
				append_cigar(aln, {{'D', anchors[*cur].query_start - anchors[*prev].query_end}});	
			} else {
				assert(anchors[*cur].ref_start - anchors[*prev].ref_end);
				append_cigar(aln, {{'I', anchors[*cur].ref_start - anchors[*prev].ref_end}});	
			}
			append_cigar(aln, alignments[*cur].cigar);

			prev = cur;
		}
		aln.calculate_error();
		hits.push_back({
			query_hash.seq, anchors[dq.front()].query_start, anchors[dq.back()].query_end,
			ref_hash.seq, anchors[dq.front()].ref_start, anchors[dq.back()].ref_end,
			9999, "", "",
			aln
		});
	}

	return hits;
}

/******************************************************************************/

auto find_chains(vector<pair<Anchor, bool>> &anchors)
{
	auto pless = [&](int i, int j) { // is anchors[i] < anchors[j] ?
		auto &pp = anchors[i].first;
		auto &p  = anchors[j].first;
		return pp.query_end < p.query_start && pp.ref_end < p.ref_start;
	};

	assert(anchors.size());

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

		int gap = p.query_start - pp.query_end;
		if (gap > MAXGAP && double(gap) / (p.query_end - pp.query_start) > MAX_ERROR) 
			return 2;
		gap = p.ref_start - pp.ref_end;
		if (gap > MAXGAP && double(gap) / (p.ref_end - pp.ref_start) > MAX_ERROR) 
			return 1;
		
		return 0;
	};
	eprn("inti S={}", S.size());
	for (int i = 1; i < path.size(); i++) {
		if (pgap(path[i - 1], path[i])) {
			splits.push_back(i);
		}
	}
	splits.push_back(path.size());

	vector<vector<int>> paths;
	for (int i = 1; i < splits.size(); i++) {
		if (anchors[path[splits[i] - 1]].first.query_start - anchors[path[splits[i - 1]]].first.query_start + KMER_SIZE > 1000 && 
			anchors[path[splits[i] - 1]].first.ref_start   - anchors[path[splits[i - 1]]].first.ref_start   + KMER_SIZE > 1000) 
		{
			paths.push_back(vector<int>());
			for (int j = splits[i - 1]; j < splits[i]; j++) {
				paths.back().push_back(path[j]);
				anchors[path[j]].second = false;
			}
		}
	}

	return paths;
}

auto cluster(const Index &query_hash, const Index &ref_hash) 
{
	auto T = cur_time();
	eprn("-- clustering len {:n} vs {:n}", query_hash.seq->seq.size(), ref_hash.seq->seq.size());

	// 1. Generate the list of hits (small anchors) inside the dot graph
	vector<pair<Anchor, bool>> anchors; // sorted by query_start, then ref_start
	for (auto &query: query_hash.minimizers) {
		auto ptr = ref_hash.index.find(query.hash);
		if (ptr == ref_hash.index.end()) 
			continue;
		for (auto ref_loc: ptr->second) {
			anchors.push_back({{
				query.loc, query.loc + query_hash.kmer_size,
				ref_loc, ref_loc + ref_hash.kmer_size, 
				{query.loc}, {ref_loc}
			}, true});
			if (anchors.size() > 1) {
				assert(anchors[anchors.size() - 2].first < anchors.back().first);
			}
		}
	}
	eprn("size: {}", anchors.size());

	// 2. Run DP on the anchors and collect all different anchors
	vector<Anchor> long_chains;
	for (int iter = 0; iter < 20 && anchors.size(); iter++) {
		auto chains = find_chains(anchors);
		if (chains.size() == 0) 
			break;
		for (auto &chain: chains) {
			int qlo = anchors[chain.front()].first.query_start, 
			    qhi = anchors[chain.back()].first.query_end;
			int rlo = anchors[chain.front()].first.ref_start,   
			    rhi = anchors[chain.back()].first.ref_end;
			Anchor a { qlo, qhi, rlo, rhi, {}, {} };
			for (int ai: chain) {
				a.query_kmers.push_back(anchors[ai].first.query_kmers.front());
				a.ref_kmers.push_back(anchors[ai].first.ref_kmers.front());
			}
			for (int ai = chain.front(); ai <= chain.back(); ai++) {
				assert(anchors[ai].first.query_start >= qlo);
				assert(anchors[ai].first.query_end <= qhi);
				if (anchors[ai].first.ref_start >= rlo && anchors[ai].first.ref_end <= rhi) {
					anchors[ai].second = false;
				}
			}
			long_chains.push_back(a);
		}

		anchors.erase(
			remove_if(anchors.begin(), anchors.end(), [](auto &x) { return !x.second; }), 
			anchors.end()
		);
		eprn(":: iter = {}, elapsed/dp = {}s", iter, elapsed(T)), T=cur_time();
	}

	eprn("-- initial chains {}", long_chains.size());
	auto hits = merge_long_chains(long_chains, query_hash, ref_hash);
	eprn("-- final chains {}", hits.size());
	eprn(":: elapsed/long = {}s", elapsed(T)), T=cur_time();
	for (auto &hit: hits) {
		eprn("||> Q {:7n}..{:7n} R {:7n}..{:7n} | L {:7n} {:7n} | E {:4.2f} (g={:4.2f}, m={:4.2f})", 
			hit.query_start, hit.query_end,
			hit.ref_start, hit.ref_end,
			hit.query_end - hit.query_start, hit.ref_end - hit.ref_start,
			hit.aln.error.error(), hit.aln.error.gap_error(), hit.aln.error.mis_error()
		);
	}
	eprn(":: elapsed/alignment = {}s", elapsed(T)), T=cur_time();
	
	sort(hits.begin(), hits.end(), [](const Hit &a, const Hit &b) {
		return max(a.query_end - a.query_start, a.ref_end - a.ref_start) > max(b.query_end - b.query_start, b.ref_end - b.ref_start);
	});
	return hits;
}

/******************************************************************************/

vector<Hit> fast_align(const string &sa, const string &sb)
{
	Index a(make_shared<Sequence>("A", sa), KMER_SIZE, WINDOW_SIZE, false);
	Index b(make_shared<Sequence>("B", sb), KMER_SIZE, WINDOW_SIZE, false);
	return cluster(a, b);
}

void test(int, char** argv)
{
	FastaReference fr("data/hg19/hg19.fa");
	// 200sec
	string s = "chr1	145291803	145450964	chr1	148576062	148687464		-1	+	+	159161	0		482	OK";
	// ifstream fin(argv[0]);
	while (1) {
		Hit h = Hit::from_bed(s);
		
		auto T = cur_time();
		const int OX = 0;
		auto q = fr.get_sequence(h.query->name, h.query_start, h.query_end);
		auto r = fr.get_sequence(h.ref->name, h.ref_start, h.ref_end);
		fast_align(q, r);
		eprn("total {}s", elapsed(T));
		break;
	}
}

// ATAGCTAGCTAGCAT
// --AGCTAcC--GCATA