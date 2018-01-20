/// 786

#include <bits/stdc++.h>

#include "common.h"
#include "fasta.h"
#include "align.h"
#include "search.h"
#include "segment.h"

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

template<typename T>
void iter(SegmentTree<T> &x, int i=0)
{
	if(i>=x.tree.size()) return;
	if (x.tree[i].a!=-1) {
		eprn("tree: {}::{}", x.anchors[x.tree[i].a].x.first, x.anchors[x.tree[i].a].x.second);
		// if(2*i+1<x.tree.size()) {
		// 	eprn("--> tree: {}::{}", x.anchors[x.tree[2*i+1].a].x.first, x.anchors[x.tree[2*i+1].a].x.second);
		// 	assert(0);
		// 	// ...
		// }
		return;
	}
	iter(x,2*i+1);
	iter(x,2*i+2);
}

auto find_chains(vector<pair<Anchor, bool>> &anchors)
{
	// QUERY is coordinate 1; REF is coordinate 2 for RMQ 
	vector<pair<int, int>> xs; // QUERY; pos in ANCHORS

	struct Coor {
		pair<int, int> x; 
		int score, pos;
		bool operator<(const Coor &a) const { return x < a.x; }
	};
	vector<Coor> ys;
	int max_q = 0, max_r = 0;
	for (int i = 0; i < anchors.size(); i++) {
		anchors[i].second = 0;
		xs.push_back({anchors[i].first.query_start, i});
		xs.push_back({anchors[i].first.query_end, i});
		ys.push_back({{anchors[i].first.ref_end - 1, i}, -1, i}); 
		max_q = max(max_q, anchors[i].first.query_end);
		max_r = max(max_r, anchors[i].first.ref_end);
	}
	sort(xs.begin(), xs.end());

	SegmentTree<Coor> tree(ys); // TODO: ys is copied; maybe change!
	// iter(tree);



	vector<int> prev(anchors.size(), -1);
	vector<int> dp(anchors.size(), 0);
	int max = 0, maxi = -1;
	for (auto &x: xs) {
		int i = x.second;
		if (x.first == anchors[i].first.query_start) {
			anchors[i].second = 1;
			int w = anchors[i].first.query_end - anchors[i].first.query_start;
			int j = tree.rmq({anchors[i].first.ref_start - 1, anchors.size()});

			if (j != -1) {
				j = ys[j].pos;
				prev[i] = j;
				int gap = anchors[i].first.query_start - anchors[j].first.query_end +
				          anchors[i].first.ref_start - anchors[j].first.ref_end;
				dp[i] = w + dp[j] - gap;
			} else {
				dp[i] = w;
			}
			if (dp[i] > max) {
				max = dp[i], maxi = i;
			}
		} else {
			assert(anchors[i].second);
			// eprn("activating {}/{}", anchors[i].first.ref_end - 1, i);
			// bool OK=0;
			// for (int W = 0; W < ys.size(); W++)
			// 	if(ys[W].x==make_pair(anchors[i].first.ref_end - 1, i)) {
			// 		OK=1; break;
			// 	}
			// assert(OK);
			int gap = max_q + 1 - anchors[i].first.query_end +
				      max_r + 1 - anchors[i].first.ref_end;
			tree.activate({anchors[i].first.ref_end - 1, i}, dp[i] - gap);
		}
	}

	assert(maxi != -1);
	deque<int> path;
	while (maxi != -1) {
		path.push_front(maxi);
		maxi = prev[maxi];
	}

	// eprn("size: {}", path.size());
	// for (int pi = 0; pi < path.size(); pi++) {
	// 	auto &a = anchors[path[pi]].first;
	// 	if (pi) {
	// 		auto &p = anchors[path[pi-1]].first;
	// 		eprn("    {} {}", a.query_start-p.query_end, a.ref_start - p.ref_end);
	// 	}
	// 	eprn("aa: {}: {} --> {}", abs(a.query_start-a.query_end), a.query_start, a.ref_start);
	// }

	return vector<vector<int>>();
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
	// for (auto &af: anchors) {
	// 	auto &a = af.first;
	// 	eprn("{}: {} --> {}", abs(a.query_start-a.query_end), a.query_start, a.ref_start);
	// }

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
	return vector<Hit>();
	// auto hits = merge_long_chains(long_chains, query_hash, ref_hash);
	// eprn("-- final chains {}", hits.size());
	// eprn(":: elapsed/long = {}s", elapsed(T)), T=cur_time();
	// for (auto &hit: hits) {
	// 	eprn("||> Q {:7n}..{:7n} R {:7n}..{:7n} | L {:7n} {:7n} | E {:4.2f} (g={:4.2f}, m={:4.2f})", 
	// 		hit.query_start, hit.query_end,
	// 		hit.ref_start, hit.ref_end,
	// 		hit.query_end - hit.query_start, hit.ref_end - hit.ref_start,
	// 		hit.aln.error.error(), hit.aln.error.gap_error(), hit.aln.error.mis_error()
	// 	);
	// }
	// eprn(":: elapsed/alignment = {}s", elapsed(T)), T=cur_time();
	
	// sort(hits.begin(), hits.end(), [](const Hit &a, const Hit &b) {
	// 	return max(a.query_end - a.query_start, a.ref_end - a.ref_start) > max(b.query_end - b.query_start, b.ref_end - b.ref_start);
	// });
	// return hits;
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