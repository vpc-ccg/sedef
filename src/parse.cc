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
};

/******************************************************************************/

auto generate_anchors(const string &query, const string &ref, const int kmer_size = 10)
{
	const uint32_t MASK = (1 << (2 * kmer_size)) - 1;

	unordered_map<uint32_t, list<int>> ref_hashes;

	int last_n = - kmer_size;
	uint32_t h = 0;
	for (int i = 0; i < ref.size(); i++) {
		if (ref[i] == 'N') 
			last_n = i;
		h = ((h << 2) | hash_dna(ref[i])) & MASK; 
		if (i < kmer_size) 
			continue;
		if (last_n >= (i - kmer_size + 1)) 
			continue;
		ref_hashes[h].push_back(i - kmer_size + 1);
	}
	
	unordered_map<int, vector<Anchor>> slide;
	last_n = -kmer_size, h = 0;
	for (int i = 0; i < query.size(); i++) {
		if (query[i] == 'N') 
			last_n = i;
		h = ((h << 2) | hash_dna(query[i])) & MASK; 
		if (i < kmer_size) 
			continue;
		if (last_n >= (i - kmer_size + 1)) 
			continue;
		
		auto it = ref_hashes.find(h);
		if (it == ref_hashes.end()) 
			continue;

		int q = i - kmer_size + 1;
		for (int r: it->second) {
			auto &diagonal = slide[r - q];
			assert(!diagonal.size() || (diagonal.back().query_start < q && diagonal.back().ref_start < r));
			if (diagonal.size() && diagonal.back().query_end >= q) {
				diagonal.back().query_end = max(diagonal.back().query_end, q + kmer_size);
				diagonal.back().ref_end = max(diagonal.back().ref_end, r + kmer_size);
			} else {
				diagonal.push_back({q, q + kmer_size, r, r + kmer_size});
			}
		}
	}
	
	vector<pair<Anchor, bool>> anchors;
	for (auto &diagonal: slide) 
		for (auto &a: diagonal.second) {
			a.query_kmers.push_back({a.query_start, a.query_end});
			a.ref_kmers.push_back({a.ref_start, a.ref_end});
			anchors.push_back({a, true});
		}
	sort(anchors.begin(), anchors.end());

	return anchors; 
}

auto chain_anchors(vector<pair<Anchor, bool>> &anchors)
{
	struct Coor {
		pair<int, int> x; 
		int score, pos;
		bool operator<(const Coor &a) const { return x < a.x; }
	};
	vector<Coor> xs, ys; // QUERY; pos in ANCHORS
	
	int max_q = 0, max_r = 0;
	for (int i = 0; i < anchors.size(); i++) {
		if (!anchors[i].second) 
			continue;
		auto &a = anchors[i].first;
		xs.push_back({{a.query_start, i}, -1, i});
		xs.push_back({{a.query_end, i}, -1, i});
		ys.push_back({{a.ref_end - 1, i}, -1, i}); 
		
		assert(a.query_start < a.query_end);
		max_q = max(max_q, a.query_end);
		max_r = max(max_r, a.ref_end);
	}

	sort(xs.begin(), xs.end());
	SegmentTree<Coor> tree(ys); 

	vector<int> prev(anchors.size(), -1);
	vector<int> dp(anchors.size(), 0);
	int max = 0, maxi = -1;

	for (auto &x: xs) {
		int i = x.x.second;
		auto &a = anchors[i].first;
		if (x.x.first == a.query_start) {
			int w = a.query_end - a.query_start;
			int j = tree.rmq({a.ref_start - 1, anchors.size()});
			if (j != -1 && ys[j].score != -1) {
				j = ys[j].pos;
				auto &p = anchors[j].first;
				assert(a.query_start >= p.query_end);
				assert(a.ref_start >= p.ref_end);
				int gap = a.query_start - p.query_end + a.ref_start - p.ref_end;
				if (w + dp[j] - gap > 0) {
					dp[i] = w + dp[j] - gap;
					prev[i] = j;
				} else {
					dp[i] = w;
				}
			} else {
				dp[i] = w;
			}
			if (dp[i] > max) {
				max = dp[i], maxi = i;
			}
		} else {
			int gap = max_q + 1 - a.query_end + max_r + 1 - a.ref_end;
			tree.activate({a.ref_end - 1, i}, dp[i] - gap);
		}
	}

	deque<int> path;
	while (maxi != -1) {
		path.push_front(maxi);
		anchors[maxi].second = false;
		maxi = prev[maxi];
	}
	// for (auto i: path) {
	// 	eprn("=== [{:4}] {:6}..{:6} -> {:6}..{:6}", i, anchors[i].first.query_start, anchors[i].first.query_end,
	// 		anchors[i].first.ref_start, anchors[i].first.ref_end);
	// }
	return path;
}

vector<Hit> fast_align(const string &query, const string &ref)
{
	auto T = cur_time();
	eprn("-- aligning query {:n} --> ref {:n}", query.size(), ref.size());

	// 1. Generate the list of hits (small anchors) inside the dot graph	
	auto anchors = generate_anchors(query, ref);

	// 2. Run DP on the anchors and collect all different anchors
	vector<Anchor> chains;
	for (int iter = 0; iter < 20; iter++) {
		auto chain = chain_anchors(anchors);
		if (chain.size() == 0) 
			break;
		
		int qlo = anchors[chain.front()].first.query_start, 
			qhi = anchors[chain.back()].first.query_end;
		int rlo = anchors[chain.front()].first.ref_start,   
			rhi = anchors[chain.back()].first.ref_end;

		if (min(rhi - rlo, qhi - qlo) < MIN_READ_SIZE)
			break;

		assert(qhi <= query.size());
		assert(rhi <= ref.size());

		Anchor a { qlo, qhi, rlo, rhi, {}, {} };
		for (int ai: chain) {
			a.query_kmers.push_back(anchors[ai].first.query_kmers.front());
			a.ref_kmers.push_back(anchors[ai].first.ref_kmers.front());
		}
		
		for (int ai = 1; ai < chain.size(); ai++)
			assert(anchors[ai - 1] < anchors[ai]);
		for (int ai = chain.front(); ai <= chain.back(); ai++) {
			assert(anchors[ai].first.query_start >= qlo);
			if (anchors[ai].first.query_end <= qhi && anchors[ai].first.ref_start >= rlo && anchors[ai].first.ref_end <= rhi) {
				anchors[ai].second = false;
			}
		}
		chains.push_back(a);
		eprn("-- chain: (len:{}) {}..{} --> {}..{}", abs(a.query_start-a.query_end), 
			a.query_start, a.query_end, a.ref_start, a.ref_end);

		eprn(":: iter = {}, elapsed/dp = {}s", iter, elapsed(T)), T=cur_time();
	}

	// 3. Perform the full alignment
	eprn("-- initial chains {}", chains.size());
	eprn(":: elapsed/long = {}s", elapsed(T)), T=cur_time();
	vector<Hit> hits;

	auto query_ptr = make_shared<Sequence>("QRY", query);
	auto ref_ptr = make_shared<Sequence>("REF", ref);
	for (auto &ch: chains) {
		auto hit = Hit {
			query_ptr, ch.query_start, ch.query_end,
			ref_ptr, ch.ref_start, ch.ref_end,
			0, "", "", {}
		};
		hit.aln = Alignment::from_anchors(query, ref, ch.query_kmers, ch.ref_kmers);
		eprn("||> Q {:7n}..{:7n} R {:7n}..{:7n} | L {:7n} {:7n} | E {:4.2f} (g={:4.2f}, m={:4.2f})", 
			hit.query_start, hit.query_end,
			hit.ref_start, hit.ref_end,
			hit.query_end - hit.query_start, hit.ref_end - hit.ref_start,
			hit.aln.error.error(), hit.aln.error.gap_error(), hit.aln.error.mis_error()
		);
		hits.push_back(hit);
	}
	eprn(":: elapsed/alignment = {}s", elapsed(T)), T=cur_time();

	return hits;
}

/******************************************************************************/

void test(int, char** argv)
{
	FastaReference fr("data/hg19/hg19.fa");
	// 200sec
	string s;
	// s = "chr22	16239131	16243489	chr22	16244049	16248400	align_both/0015/both076775	-1.0	+	+	4358	0		-nan	err=7.5;mis=-nan;gap=-nan";
	s = "chr1	145883118	146164650	chr1	147424817	147706477	NA	NA	+	+";
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