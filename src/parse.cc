/// 786

#include <bits/stdc++.h>

#include "common.h"
#include "fasta.h"
#include "align.h"
#include "search.h"

using namespace std;

/******************************************************************************/

// #define eprn

struct Anchor {
	int query_start, query_end;
	int ref_start, ref_end;
	list<pair<int, int>> query_kmers, ref_kmers;

	bool operator< (const Anchor &a) const {
		return tie(query_start, ref_start) < tie(a.query_start, a.ref_start);
	}
	Alignment anchor_align(const string &qs, const string &rs) ;
};

template<typename T>
struct SegmentTree {
	struct Point {
		int p, h, a; /* h is inclusive; a is index in the original array; -1 otherwise */
	};
	vector<Point> tree;
	vector<T> anchors;

	int rmq(int q, int i, int s, int e) // RMQ for [0, q) in [s, e): returns index i in tree
	{
		if (i >= tree.size()) {
			return -1;
		} else if (s + 1 == e) { // leaf
			if (anchors[tree[i].a].x < q)
				return i;
		} else {
			assert(tree[i].p != -1);
			int &pv = tree[i].p;
			assert(tree[pv].a != -1);
			if (anchors[tree[pv].a].x < q) {
				return pv;
			} else {
				assert(2 * i + 1 < tree.size());
				if (q - 1 <= tree[2 * i + 1].h) { // h is inclusive
					return rmq(q, 2 * i + 1, s, (s + e) / 2);
				} else {
					int m1 = rmq(q, 2 * i + 1, s, (s + e) / 2);
					int m2 = rmq(q, 2 * i + 2, (s + e) / 2, e); 
					if (m1 == -1) return m2;
					if (m2 == -1) return m1;
					return anchors[tree[m1].a].score >= anchors[tree[m2].a].score ? m1 : m2;
				}
			}
		}
		return -1;
	}
	T& rmq(int q)
	{
		int i = rmq(q, 0, 0, tree.size());
		assert(i != -1);
		return anchors[tree[i].a];
	}

	void activate(int leaf, int i, int s, int e)
	{
		if (tree[i].p == -1 || anchors[tree[tree[i].p].a].score <= anchors[tree[leaf].a].score)
			swap(tree[i].p, leaf);
		if (leaf == -1) 
			return;

		assert(2 * i + 1 < tree.size());
		if (anchors[tree[leaf].a].x <= tree[2 * i + 1].h) {
			activate(leaf, 2 * i + 1, s, (s + e) / 2);
		} else {
			assert(2 * i + 2 < tree.size());
			activate(leaf, 2 * i + 2, (s + e) / 2, e);
		}
	}

	void activate(int q)
	{
		// eprn("looking for {}", q);
		int i;
		for (i = 0; i < tree.size() && q != anchors[tree[i].a].x; ) {
			i = 2 * i + 1 + (q > tree[2 * i + 1].h);
		}
		assert(i < tree.size());
		eprn("found {} at {}", q, i);
		activate(i, 0, 0, tree.size());
	}

	int tree_i;
	int initialize(int i, int s, int e)
	{
		assert(i < tree.size());
		if (s + 1 == e) {
			assert(tree_i < anchors.size());
			tree[i] = Point {-1, anchors[tree_i].x, tree_i};
			tree_i++;
			return i;
		} else {
			int a = initialize(2 * i + 1, s, (s + e) / 2);
			int b = initialize(2 * i + 2, (s + e) / 2, e);
			tree[i] = {-1, 2 * i + 2 < tree.size() ? tree[2 * i + 2].h : tree[2 * i + 1].h, -1};
			return max(a, b);
		}
	}

	SegmentTree(const vector<T> &a): anchors(a)
	{
		sort(anchors.begin(), anchors.end());
		tree.resize(3 * anchors.size());

		tree_i = 0;
		int m = initialize(0, 0, anchors.size());
		while (tree.size() > m + 1) tree.pop_back();
		assert(tree.back().a == anchors.size() - 1);
	}


	void plot(int l, int i, int s, int e, vector<string> &PLOT)
	{
		if (i>=tree.size()) return;
		plot(l+1,2*i+1,s,(s+e)/2,PLOT);
		PLOT[l]+=fmt::format("{}::{} ",
			tree[i].a==-1?-1:anchors[tree[i].a].x, 
			tree[i].p != -1 ? &anchors[tree[tree[i].p].a]-&anchors[0]+1 : -1);
		plot(l+1,2*i+2,(s+e)/2,e,PLOT);
	}
	void plot()  
	{
		vector<string> PLOT(50);
		plot(0, 0, 0, tree.size(),PLOT);
		for (auto &s: PLOT) if (s=="") break; else eprn("{}: {}", &s-&PLOT[0], s);
	}
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

Alignment Anchor::anchor_align(const string &qstr, const string &rstr) 
{
	auto aln_from_match = [&](int qs, int qe, int rs, int re) {
		assert(qe - qs == re - rs);
		return Alignment {
			"A", qs, qe, 
			"B", rs, re,
			qstr.substr(qs, qe - qs),
			rstr.substr(rs, re - rs),
			"", "", "",
			{{'M', qe - qs}}, {}
		};
	};

	auto qk_p = query_kmers.begin();
	auto rk_p = ref_kmers.begin();
	Alignment aln = aln_from_match(qk_p->first, qk_p->second, rk_p->first, rk_p->second);
	assert(query_kmers.size() == ref_kmers.size());
	for (auto qk = next(qk_p), rk = next(rk_p); qk != query_kmers.end(); qk++, rk++) {
		eprn("--- {}..{} vs {}..{}", qk_p->first, qk_p->second,
			rk_p->first, rk_p->second);
		aln.end_a = qk->second;
		aln.end_b = rk->second;
		aln.a += qstr.substr(qk_p->second, qk->second - qk_p->second);
		aln.b += rstr.substr(rk_p->second, rk->second - rk_p->second);
		if (qk->first - qk_p->second && rk->first - rk_p->second) {
			auto gap = align(
				qstr.substr(qk_p->second, qk->first - qk_p->second), 
				rstr.substr(rk_p->second, rk->first - rk_p->second), 
				5, -4, 40, 1
			);
			append_cigar(aln, gap.cigar);
		} else if (qk->first - qk_p->second) {
			append_cigar(aln, {{'D', qk->first - qk_p->second}});	
		} else if (rk->first - rk_p->second) {
			append_cigar(aln, {{'I', rk->first - rk_p->second}});	
		} 
		assert(qk->second - qk->first == rk->second - rk->first);
		append_cigar(aln, {{'M', qk->second - qk->first}});
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
		auto &aln = alignments.back();
		// eprn("{}..{} -> {}..{}", aln.start_a, aln.end_a, aln.start_b, aln.end_b);
		// eprn("{}", aln.print());
	}
	exit(0);

	vector<int> score, dp; // fill the scores!!
	for (int i = 0; i < anchors.size(); i++) {
		auto &err = alignments[i].error;
		score.push_back(err.matches * 5 - err.mismatches * 4 - 40 * err.gaps - err.gap_bases);
		dp.push_back(score[i]);
	}
	vector<int> prev(anchors.size(), -1);
	// maximize span while keeping the properties
	// for (int i = 0; i < anchors.size(); i++) {
	// 	for (int j = 0; j < i; j++) {
	// 		if (anchors[j].query_end >= anchors[i].query_start) 
	// 			continue;
	// 		if (anchors[j].ref_end >= anchors[i].ref_start) 
	// 			continue;

	// 		int h = anchors[i].query_start - anchors[j].query_end;
	// 		int v = anchors[i].ref_start - anchors[j].ref_end;
	// 		int dist = score[i] - (4 * min(h, v) + 40 + abs(h - v));

	// 		if (dp[j] + dist > dp[i]) {
	// 			dp[i] = dp[j] + dist;
	// 			prev[i] = j;
	// 		}
	// 	}
	// }

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
		return pp.query_end <= p.query_start && pp.ref_end <= p.ref_start;
	};

	assert(anchors.size());

	vector<int> prev(anchors.size(), -1);
	vector<int> S { 0 };
	for (int i = 1; i < anchors.size(); i++) {
		if (pless(S.back(), i)) {
			prev[i] = S.back();
			S.push_back(i);
		} else {
			// returns first i s.t. X[i] >= y
			auto pos = distance(S.begin(), lower_bound(S.begin(), S.end(), i, [&](int i, int j) {
				return anchors[i].first < anchors[j].first;
			}));
			assert(pos >= 0 && pos < S.size());
			prev[i] = pos ? S[pos - 1] : -1; 
			S[pos] = i;
		}
	}

	deque<int> path;
	vector<int> splits{0};

	int cur = S.back();
	for (int i = 0; i < S.size(); i++) {
		assert(cur != -1);
		path.push_front(cur);
		cur = prev[cur];
	}

	eprn("size: {}", S.size());
	for (auto &p: path) {
		auto &a = anchors[p].first;
		eprn("aa: {}: {} --> {}", abs(a.query_start-a.query_end), a.query_start, a.ref_start);
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
	for (int i = 1; i < path.size(); i++) {
		if (pgap(path[i - 1], path[i])) {
			// splits.push_back(i);
		}
	}
	splits.push_back(path.size());

	vector<vector<int>> paths;
	for (int i = 1; i < splits.size(); i++) {
		if (anchors[path[splits[i] - 1]].first.query_start - anchors[path[splits[i - 1]]].first.query_end > 1000 && 
			anchors[path[splits[i] - 1]].first.ref_start   - anchors[path[splits[i - 1]]].first.query_end > 1000) 
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


auto chain_full(vector<Anchor> &anchors) 
{
	vector<pair<int, int>> xes; // query_s, q_e

	for (int i = 0; i < anchors.size(); i++) {
		xes.push_back({an.query_start, i});
		xes.push_back({an.query_end, i});	
	}
	sort(xes.begin(), xes.end());

	auto tree = segment_tree(anchors);
	for (auto &x: xes) { // querys
		int i = x.second;
		if (x.first == anchors[i].query_start) {
			int q = tree.rmq(anchors[i].ref_start, 0, 0, tree.size());
			int j = tree.tree[q].anchor;
			prev[i] = j;
			dp[i] = weight(anchors[i]) + dp[j];
		} else {
			activate("ANCHOOOORRRRRR", 0, 0, ) ;
		}
	}
}

auto cluster(const Index &query_hash, const Index &ref_hash) 
{
	auto T = cur_time();
	eprn("-- clustering len {:n} vs {:n}", query_hash.seq->seq.size(), ref_hash.seq->seq.size());

	// 1. Generate the list of hits (small anchors) inside the dot graph
	int kmer_size = 10;
	const uint32_t MASK = (1 << (2 * kmer_size)) - 1;

	unordered_map<uint32_t, list<int>> ref;
	int last_n = - kmer_size;
	uint32_t h = 0;
	for (int i = 0; i < ref_hash.seq->seq.size(); i++) {
		if (ref_hash.seq->seq[i] == 'N') last_n = i;
		h = ((h << 2) | hash_dna(ref_hash.seq->seq[i])) & MASK; 
		if (i < kmer_size) continue;
		if (last_n >= (i - kmer_size + 1)) continue;
		ref[h].push_back(i);
	}
	eprn("ref hash {}", ref.size());
	unordered_map<int, vector<Anchor>> slide;
	int ss = 0, sq = 0;
	last_n = -kmer_size;
	h = 0;
	for (int i = 0; i < query_hash.seq->seq.size(); i++) {
		if (query_hash.seq->seq[i] == 'N') last_n = i;
		h = ((h << 2) | hash_dna(query_hash.seq->seq[i])) & MASK; 
		if (i < kmer_size) continue;
		if (last_n >= (i - kmer_size + 1)) continue;
		
		auto it = ref.find(h);
		if (it == ref.end()) continue;

		int q = i;
		for (int r: it->second) {
			auto &lii = slide[r - q];
			ss++;
			if (lii.size()) {
				assert(lii.back().query_start < q && lii.back().ref_start < r);
			}
			if (lii.size() && lii.back().query_end >= q) {
				lii.back().query_end = max(lii.back().query_end, q + kmer_size);
				lii.back().ref_end = max(lii.back().ref_end, r + kmer_size);
			} else {
				sq++;
				lii.push_back({q, q + kmer_size, r, r + kmer_size});
			}
		}
	}
	eprn("init {}, merge {}", ss, sq);

	vector<pair<Anchor, bool>> anchors; // sorted by query_start, then ref_start
	for (auto &sx: slide) 
		for (auto &ao: sx.second) {

			anchors.push_back({ao, true});
			anchors.back().first.query_kmers.push_back({anchors.back().first.query_start, anchors.back().first.query_end});
			anchors.back().first.ref_kmers.push_back({anchors.back().first.ref_start, anchors.back().first.ref_end});
		}
	sort(anchors.begin(), anchors.end());

	eprn("size: {}", anchors.size());
	for (auto &af: anchors) {
		auto &a = af.first;
		eprn("{}: {} --> {}", abs(a.query_start-a.query_end), a.query_start, a.ref_start);
	}

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
			eprn("aa: {}: {} --> {}", abs(a.query_start-a.query_end), a.query_start, a.ref_start);

		}
		exit(0);

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

// #include <seqan/seeds.h>
// auto seqqan(const Index &query_hash, const Index &ref_hash) 
// {
// 	using namespace seqan;

// 	auto T = cur_time();
// 	eprn("-- clustering len {:n} vs {:n}", query_hash.seq->seq.size(), ref_hash.seq->seq.size());

// 	// 1. Generate the list of hits (small anchors) inside the dot graph
// 	int ss = 0;
// 	SeedSet<Seed<Simple>> anchors; 
// 	for (auto &query: query_hash.minimizers) {
// 		auto ptr = ref_hash.index.find(query.hash);
// 		if (ptr == ref_hash.index.end()) 
// 			continue;
// 		for (auto ref_loc: ptr->second) {
// 			auto seed = Seed<Simple>(
// 				query.loc, ref_loc,
// 				query.loc + query_hash.kmer_size,
// 				ref_loc + ref_hash.kmer_size);
// 			ss++;
// 			if (!addSeed(anchors, seed, 0, 0, Score<int, Simple>(), 0, 0, Merge()))
// 				addSeed(anchors, seed, Single());
// 		}
// 	}
// 	eprn("size: {} ({})", length(anchors), ss);
// 	eprn(":: elapsed/building = {}s", elapsed(T)), T=cur_time();


// 	String<Seed<Simple>> result;
// 	chainSeedsGlobally(result, anchors, SparseChaining());
// 	eprn(":: elapsed/chaining = {}s", elapsed(T)), T=cur_time();
	
// 	Align<Dna5String, ArrayGaps> alignment;
// 	resize(rows(alignment), 2);
// 	assignSource(row(alignment, 0), query_hash.seq->seq);
// 	assignSource(row(alignment, 1), ref_hash.seq->seq);

// 	Score<int, Simple> scoring(5, -4, -1, -40);
// 	int r = bandedChainAlignment(alignment, result, scoring);
// 	eprn(":: elapsed/alignment = {}s", elapsed(T)), T=cur_time();

// 	// eprn("score={}\n{}", r, alignment);

// 	// for (auto &af: anchors) {
// 	// 	auto &a = af.first;
// 	// 	eprn("{}: {} --> {}", abs(a.query_start-a.ref_start), a.query_start, a.ref_start);
// 	// }

// 	// 2. Run DP on the anchors and collect all different anchors
// 	// vector<Anchor> long_chains;
// 	// for (int iter = 0; iter < 20 && anchors.size(); iter++) {
// 	// 	auto chains = find_chains(anchors);
// 	// 	if (chains.size() == 0) 
// 	// 		break;
// 	// 	for (auto &chain: chains) {
// 	// 		int qlo = anchors[chain.front()].first.query_start, 
// 	// 		    qhi = anchors[chain.back()].first.query_end;
// 	// 		int rlo = anchors[chain.front()].first.ref_start,   
// 	// 		    rhi = anchors[chain.back()].first.ref_end;
// 	// 		Anchor a { qlo, qhi, rlo, rhi, {}, {} };
// 	// 		for (int ai: chain) {
// 	// 			a.query_kmers.push_back(anchors[ai].first.query_kmers.front());
// 	// 			a.ref_kmers.push_back(anchors[ai].first.ref_kmers.front());
// 	// 		}
// 	// 		for (int ai = chain.front(); ai <= chain.back(); ai++) {
// 	// 			assert(anchors[ai].first.query_start >= qlo);
// 	// 			assert(anchors[ai].first.query_end <= qhi);
// 	// 			if (anchors[ai].first.ref_start >= rlo && anchors[ai].first.ref_end <= rhi) {
// 	// 				anchors[ai].second = false;
// 	// 			}
// 	// 		}
// 	// 		long_chains.push_back(a);
// 	// 	}

// 	// 	anchors.erase(
// 	// 		remove_if(anchors.begin(), anchors.end(), [](auto &x) { return !x.second; }), 
// 	// 		anchors.end()
// 	// 	);
// 	// 	eprn(":: iter = {}, elapsed/dp = {}s", iter, elapsed(T)), T=cur_time();
// 	// }

// 	// eprn("-- initial chains {}", long_chains.size());
// 	// auto hits = merge_long_chains(long_chains, query_hash, ref_hash);
// 	// eprn("-- final chains {}", hits.size());
// 	// eprn(":: elapsed/long = {}s", elapsed(T)), T=cur_time();
// 	// for (auto &hit: hits) {
// 	// 	eprn("||> Q {:7n}..{:7n} R {:7n}..{:7n} | L {:7n} {:7n} | E {:4.2f} (g={:4.2f}, m={:4.2f})", 
// 	// 		hit.query_start, hit.query_end,
// 	// 		hit.ref_start, hit.ref_end,
// 	// 		hit.query_end - hit.query_start, hit.ref_end - hit.ref_start,
// 	// 		hit.aln.error.error(), hit.aln.error.gap_error(), hit.aln.error.mis_error()
// 	// 	);
// 	// }
// 	// eprn(":: elapsed/alignment = {}s", elapsed(T)), T=cur_time();
	
// 	// sort(hits.begin(), hits.end(), [](const Hit &a, const Hit &b) {
// 	// 	return max(a.query_end - a.query_start, a.ref_end - a.ref_start) > max(b.query_end - b.query_start, b.ref_end - b.ref_start);
// 	// });
// 	return vector<Hit>();
// }

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
	string s;
	s = "chr22	16239131	16243489	chr22	16244049	16248400	align_both/0015/both076775	-1.0	+	+	4358	0		-nan	err=7.5;mis=-nan;gap=-nan";
	// s = "chr1	145883118	146164650	chr1	147424817	147706477	NA	NA	+	+";
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