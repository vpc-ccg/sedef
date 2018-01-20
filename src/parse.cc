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
	auto qk_p = query_kmers.begin();
	auto rk_p = ref_kmers.begin();
	Alignment aln {
		"A", qk_p->first, qk_p->second, 
		"B", rk_p->first, rk_p->second,
		qstr.substr(qk_p->first, qk_p->second - qk_p->first),
		rstr.substr(rk_p->first, rk_p->second - rk_p->first),
		"", "", "",
		{{'M', qk_p->second - qk_p->first}}, {}
	};
	assert(query_kmers.size() == ref_kmers.size());

	for (auto qk = next(qk_p), rk = next(rk_p); qk != query_kmers.end(); qk++, rk++) {
		// eprnn("--- {:6}..{:6} vs {:6}..{:6}", qk_p->first, qk_p->second,
		// 	rk_p->first, rk_p->second);
		assert(qk->first >= qk_p->second);
		assert(rk->first >= rk_p->second);
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
			// eprnn(" {}", gap.cigar_string());
			append_cigar(aln, gap.cigar);
		} else if (qk->first - qk_p->second) {
			append_cigar(aln, {{'D', qk->first - qk_p->second}});	
		} else if (rk->first - rk_p->second) {
			append_cigar(aln, {{'I', rk->first - rk_p->second}});	
		} 
		// eprn("");
		assert(qk->second - qk->first == rk->second - rk->first);
		append_cigar(aln, {{'M', qk->second - qk->first}});
		qk_p = qk, rk_p = rk;
	}

	int qlo = query_kmers.front().first, 
		qhi = query_kmers.back().second;
	int rlo = ref_kmers.front().first,   
		rhi = ref_kmers.back().second;
	assert(qhi <= qstr.size());
	assert(rhi <= rstr.size());
	assert(aln.a == qstr.substr(qlo, qhi - qlo));
	assert(aln.b == rstr.substr(rlo, rhi - rlo));

	aln.populate_nice_alignment();
	aln.error = aln.calculate_error();
	return aln;
}

/******************************************************************************/

// #include <seqan/seeds.h>

// 1..281532 --> 1..281660
// 281532M1I127M

// template<typename TAlign>
// string seqan_cigar(const TAlign &align) 
// { 
//     auto & row0 = row(align, 0);
//     auto & row1 = row(align, 1);

//     string cigar;

//     int pos = 0;
//     SEQAN_ASSERT_EQ(length(row0), length(row1));
//     int dbEndPos = length(row0);
//     int queryEndPos = length(row1);

//     int readBasePos = pos + clippedBeginPosition(row1);
//     int readPos = 0;
// 	while (pos < dbEndPos || pos < queryEndPos) {
// 		int matched = 0;
// 		int inserted = 0;
// 		int deleted = 0;
// 		while (pos != dbEndPos && pos != queryEndPos && !isGap(row0, pos) && !isGap(row1, pos)) {
//             ++readPos;
// 			++readBasePos;
// 			++pos;
// 			++matched;
// 		}
// 		if (matched > 0) 
// 			cigar += fmt::format("{}M", matched);
// 		while (pos < dbEndPos && isGap(row1, pos)) {
// 			++pos;
// 			++deleted;
// 		}
// 		if (deleted > 0) 
// 			cigar += fmt::format("{}D", deleted);
// 		while (pos < queryEndPos && isGap(row0, pos)) {
// 			++pos;
// 			++readPos;
// 			++readBasePos;
// 			++inserted;
// 		}
// 		if (inserted > 0) 
// 			cigar += fmt::format("{}I", inserted);
// 	}

// 	return cigar;
// }


// auto seqan_align(const Anchor &a, const string &q, const string &r)
// {
// 	using namespace seqan;

// 	String<Seed<Simple>> anchors; 
// 	for (auto qi = a.query_kmers.begin(), ri = a.ref_kmers.begin(); qi != a.query_kmers.end(); qi++, ri++) {
// 		appendValue(anchors, Seed<Simple>(
// 			qi->first, ri->first,
// 			qi->second - qi->first
// 		));
// 	}

// 	Align<Dna5String, ArrayGaps> alignment;
// 	resize(rows(alignment), 2);
// 	assignSource(row(alignment, 0), q);
// 	assignSource(row(alignment, 1), r);

// 	Score<int, Simple> scoring(5, -4, -1, -40);
// 	String<Seed<Simple>> result;
// 	int rc = bandedChainAlignment(alignment, result, scoring);

// 	eprn("{}", alignment);

// 	string cigar = seqan_cigar(alignment);
// 	// eprn("{}", cigar);

// 	auto aln = ::Alignment::from_cigar(q, r, cigar);
// 	return aln;
// }

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
		if (!anchors[i].second) 
			continue;
		auto &a = anchors[i].first;
		xs.push_back({a.query_start, i});
		xs.push_back({a.query_end, i});
		assert(a.query_start < a.query_end);
		ys.push_back({{a.ref_end - 1, i}, -1, i}); 
		max_q = max(max_q, a.query_end);
		max_r = max(max_r, a.ref_end);
	}
	sort(xs.begin(), xs.end());

	SegmentTree<Coor> tree(ys); // TODO: ys is copied; maybe change!

	vector<int> prev(anchors.size(), -1);
	vector<int> dp(anchors.size(), 0);
	int max = 0, maxi = -1;

	// set<int> added;
	for (auto &x: xs) {
		int i = x.second;
		if (x.first == anchors[i].first.query_start) {
			int w = anchors[i].first.query_end - anchors[i].first.query_start;
			int j = tree.rmq({anchors[i].first.ref_start - 1, anchors.size()});
			if (j != -1 && ys[j].score != -1) { // not added to the tree
				j = ys[j].pos;
				// assert(added.find(j) != added.end());
				// eprn("i: {} ({}:{}) connects to j: {} ({}:{})",
				// 	i,
				// 	anchors[i].first.query_start, anchors[i].first.query_end,
				// 	j,
				// 	anchors[j].first.query_start, anchors[j].first.query_end
				// 	);
				assert(anchors[i].first.query_start >= anchors[j].first.query_end);
				assert(anchors[i].first.ref_start >= anchors[j].first.ref_end);
				int gap = anchors[i].first.query_start - anchors[j].first.query_end +
				          anchors[i].first.ref_start - anchors[j].first.ref_end;
				dp[i] = w + dp[j] - gap;
				prev[i] = j;
			} else {
				dp[i] = w;
			}
			if (dp[i] > max) {
				max = dp[i], maxi = i;
			}
		} else {
			int gap = max_q + 1 - anchors[i].first.query_end +
				      max_r + 1 - anchors[i].first.ref_end;
			tree.activate({anchors[i].first.ref_end - 1, i}, dp[i] - gap);
			// added.insert(i);
		}
	}

	deque<int> path;
	while (maxi != -1) {
		path.push_front(maxi);
		anchors[maxi].second = false;
		maxi = prev[maxi];
	}
	// for (auto i: path) {
	// 	eprn("--> {:6}..{:6} vs {:6}..{:6}", anchors[i].first.query_start, anchors[i].first.query_end,
	// 		anchors[i].first.ref_start, anchors[i].first.ref_end);
	// }
	// exit(0);
	return path;
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
		ref[h].push_back(i - kmer_size + 1);
	}
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

		int q = i - kmer_size + 1;
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
	eprn("-- init {}, merge {}", ss, sq);

	vector<pair<Anchor, bool>> anchors;
	for (auto &sx: slide) 
		for (auto &ao: sx.second) {
			anchors.push_back({ao, true});
			anchors.back().first.query_kmers.push_back({anchors.back().first.query_start, anchors.back().first.query_end});
			anchors.back().first.ref_kmers.push_back({anchors.back().first.ref_start, anchors.back().first.ref_end});
		}
	sort(anchors.begin(), anchors.end());

	// eprn("size: {}", anchors.size());
	// for (auto &af: anchors) {
	// 	auto &a = af.first;
	// 	eprn("{}: {} --> {}", abs(a.query_start-a.query_end), a.query_start, a.ref_start);
	// }

	// 2. Run DP on the anchors and collect all different anchors
	vector<Anchor> long_chains;
	for (int iter = 0; iter < 20 && anchors.size(); iter++) {
		auto chain = find_chains(anchors);
		if (chain.size() == 0) 
			break;
		
		int qlo = anchors[chain.front()].first.query_start, 
			qhi = anchors[chain.back()].first.query_end;
		int rlo = anchors[chain.front()].first.ref_start,   
			rhi = anchors[chain.back()].first.ref_end;

		if (min(rhi - rlo, qhi - qlo) < MIN_READ_SIZE / 2)
			break;

		assert(qhi <= query_hash.seq->seq.size());
		assert(rhi <= ref_hash.seq->seq.size());

		Anchor a { qlo, qhi, rlo, rhi, {}, {} };
		for (int ai: chain) {
			a.query_kmers.push_back(anchors[ai].first.query_kmers.front());
			a.ref_kmers.push_back(anchors[ai].first.ref_kmers.front());
		}
		for (int ai = 1; ai < chain.size(); ai++)
			assert(anchors[ai-1]<anchors[ai]);

		for (int ai = chain.front(); ai <= chain.back(); ai++) {
			assert(anchors[ai].first.query_start >= qlo);
			assert(anchors[ai].first.query_end <= qhi);
			if (anchors[ai].first.ref_start >= rlo && anchors[ai].first.ref_end <= rhi) {
				anchors[ai].second = false;
			}
		}
		long_chains.push_back(a);
		eprn("-- chain: ({}) {}..{} --> {}..{}", abs(a.query_start-a.query_end), 
			a.query_start, a.query_end, a.ref_start, a.ref_end);

		
		eprn(":: iter = {}, elapsed/dp = {}s", iter, elapsed(T)), T=cur_time();
	}

	eprn("-- initial chains {}", long_chains.size());
	eprn(":: elapsed/long = {}s", elapsed(T)), T=cur_time();
	vector<Hit> hits;
	for (auto &ch: long_chains) {
		auto hit = Hit {
			query_hash.seq, ch.query_start, ch.query_end,
			ref_hash.seq, ch.ref_start, ch.ref_end,
			0, "", "", {}
		};
		hit.aln = ch.anchor_align(query_hash.seq->seq, ref_hash.seq->seq);
		eprn("||> Q {:7n}..{:7n} R {:7n}..{:7n} | L {:7n} {:7n} | E {:4.2f} (g={:4.2f}, m={:4.2f})", 
			hit.query_start, hit.query_end,
			hit.ref_start, hit.ref_end,
			hit.query_end - hit.query_start, hit.ref_end - hit.ref_start,
			hit.aln.error.error(), hit.aln.error.gap_error(), hit.aln.error.mis_error()
		);
		hits.push_back(hit);
	}
	eprn(":: elapsed/alignment = {}s", elapsed(T)), T=cur_time();

	// for (auto &hit: long_chains) {
	// 	auto aln = seqan_align(hit, query_hash.seq->seq, ref_hash.seq->seq);
	// 	eprn("||> Q {:7n}..{:7n} R {:7n}..{:7n} | L {:7n} {:7n} | E {:4.2f} (g={:4.2f}, m={:4.2f})", 
	// 		hit.query_start, hit.query_end,
	// 		hit.ref_start, hit.ref_end,
	// 		hit.query_end - hit.query_start, hit.ref_end - hit.ref_start,
	// 		aln.error.error(), aln.error.gap_error(), aln.error.mis_error()
	// 	);
	// }
	// eprn(":: elapsed/seqan = {}s", elapsed(T)), T=cur_time();
	
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
	auto h = cluster(a, b);
	return h;
}

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