/// 786

#include "common.h"
#include "fasta.h"
#include "align.h"
#include "search.h"
#include "segment.h"
#include "chain.h"

using namespace std; 

/******************************************************************************/

auto generate_anchors(const string &query, const string &ref, const int kmer_size = 12)
{
	const uint32_t MASK = (1 << (2 * kmer_size)) - 1;

	unordered_map<uint32_t, list<int>> ref_hashes;

	int last_n = - kmer_size;
	uint32_t h = 0;
	for (int i = 0; i < ref.size(); i++) {
		if (ref[i] == 'N') 
			last_n = i;
		h = ((h << 2) | hash_dna(ref[i])) & MASK; 
		if (i < kmer_size - 1) 
			continue;
		if (last_n >= (i - kmer_size + 1)) 
			continue;
		ref_hashes[h].push_back(i - kmer_size + 1);
	}
	
	vector<int> slide(query.size() + ref.size(), -1);
	vector<pair<Anchor, bool>> anchors;

	last_n = -kmer_size, h = 0;
	for (int i = 0; i < query.size(); i++) {
		if (query[i] == 'N') 
			last_n = i;
		h = ((h << 2) | hash_dna(query[i])) & MASK; 
		if (i < kmer_size - 1) 
			continue;
		if (last_n >= (i - kmer_size + 1)) 
			continue;
		
		auto it = ref_hashes.find(h);
		if (it == ref_hashes.end()) 
			continue;

		int q = i - kmer_size + 1;
		int off = query.size();
		for (int r: it->second) {
			int d = off + r - q;
			assert(d >= 0 && d < slide.size());
			if (q >= slide[d]) {
				assert(r >= d - off + slide[d]);
	
				bool has_u = 0;
				int len;
				for (len = 0; q + len < query.size() && r + len < ref.size(); len++) {
					assert(query[q + len] != 'N' && ref[r + len] != 'N');
					has_u |= (isupper(query[q + len]) || isupper(ref[r + len]));
					if (toupper(query[q + len]) != toupper(ref[r + len]))
						break;
				}
				if (len >= kmer_size / 2) {
					anchors.push_back(make_pair(
						Anchor{q, q + len, r, r + len, {{q, q + len}}, {{r, r + len}}, has_u}, true));
					slide[d] = q + len;
				}
			} else {
				assert(slide[d] >= q + kmer_size);
				assert(d - off + slide[d] >= r + kmer_size);
			}
		}
	}
	
	for (int i = 1; i < anchors.size(); i++) {
		assert(anchors[i - 1] < anchors[i]);
	}
	// sort(anchors.begin(), anchors.end());

	return anchors; 
}

string XX,YY;

auto chain_anchors(vector<pair<Anchor, bool>> &anchors)
{
	struct Coor {
		pair<int, int> x; 
		int score, pos;
		bool operator<(const Coor &a) const { return x < a.x; }
	};
	vector<Coor> xs, ys; // QUERY; pos in ANCHORS
	
	int max_q = 0, max_r = 0, l = 0;
	for (int i = 0; i < anchors.size(); i++) {
		if (!anchors[i].second) 
			continue;
		l++;
		auto &a = anchors[i].first;
		xs.push_back({{a.query_start, i}, numeric_limits<int>::min(), i});
		xs.push_back({{a.query_end, i}, numeric_limits<int>::min(), i});
		ys.push_back({{a.ref_end - 1, i}, numeric_limits<int>::min(), i}); 
		
		assert(a.query_start < a.query_end);
		max_q = max(max_q, a.query_end);
		max_r = max(max_r, a.ref_end);
	}

	// for (auto a: anchors) {
	// 	dprn("!== {:6}..{:6} -> {:6}..{:6}", a.first.query_start, a.first.query_end,
	// 		a.first.ref_start, a.first.ref_end);
	// }

	dprn("-- anchors to dp: {}", l);

	sort(xs.begin(), xs.end());
	SegmentTree<Coor> tree(ys); 

	vector<int> prev(anchors.size(), -1);
	vector<int> dp(anchors.size(), 0);
	vector<int> size_so_far(anchors.size(), 0);

	set<pair<int, int>, greater<pair<int, int>>> maxes;
	const double ratio = 4;
	const int max_gap = MAX_ERROR * MIN_READ_SIZE;
	int deactivate_bound = 0;

	for (auto &x: xs) {
		int i = x.x.second;
		auto &a = anchors[i].first;
		if (x.x.first == a.query_start) {
			while (deactivate_bound < (&x - &xs[0])) {
				int t = xs[deactivate_bound].x.second; // index
				if (xs[deactivate_bound].x.first == anchors[t].first.query_end) { // end point
					if (a.query_start - anchors[t].first.query_end <= max_gap)
						break;
					tree.deactivate({anchors[t].first.ref_end - 1, t});
				}
				deactivate_bound++;
			}

			int w = ratio * (a.query_end - a.query_start);
			int j = tree.rmq({a.ref_start - max_gap, 0}, 
				             {a.ref_start - 1, anchors.size()});
			if (j != -1 && ys[j].score != numeric_limits<int>::min()) {
				j = ys[j].pos;
				auto &p = anchors[j].first;
				assert(a.query_start >= p.query_end);
				assert(a.ref_start >= p.ref_end);
				int gap = (a.query_start - p.query_end + a.ref_start - p.ref_end);
				if (w + dp[j] - gap > 0) {
					dp[i] = w + dp[j] - gap;
					prev[i] = j;
					size_so_far[i] = size_so_far[j] + (a.query_end - a.query_start) + 
						max(a.query_start - p.query_end, a.ref_start - p.ref_end);
				} else {
					dp[i] = w;
					size_so_far[i] = (a.query_end - a.query_start);
				}
			} else {
				dp[i] = w;
				size_so_far[i] = (a.query_end - a.query_start);
			}

			// if (e.query_end <=)
			// eprn("{}..{} -> {} and {} >= {}",
				// a.query_start, a.query_end, size_so_far[i], dp[i], max(MIN_READ_SIZE, size_so_far[i]) * (ratio * (1 - MAX_GAP_ERROR) - MAX_GAP_ERROR));
			if (dp[i] >= MIN_READ_SIZE * (ratio * (1 - MAX_ERROR) - MAX_ERROR)) {
				maxes.insert({dp[i], i});
			}
		} else {
			int gap = (max_q + 1 - a.query_end + max_r + 1 - a.ref_end);
			tree.activate({a.ref_end - 1, i}, dp[i] - gap);
		}
	}


	// dp=7083 exp_size=7336

	// for (auto &m: maxes) {
	// 	eprn("{}..{} [+] {} -> {} -> {} [{}]", 
	// 		anchors[m.second].first.query_start, 
	// 		anchors[m.second].first.ref_start, 
	// 		anchors[m.second].first.query_end-anchors[m.second].first.query_start,
	// 		m.first,
	// 		XX.substr(anchors[m.second].first.query_start,anchors[m.second].first.query_end-anchors[m.second].first.query_start),
	// 		m.second
	// 	);
	// }

	vector<deque<int>> paths;
	for (auto &m: maxes) {
		int maxi = m.second;
		if (!anchors[maxi].second)
			continue;
		paths.push_back(deque<int>());


		// bool foooonf = 0;
		while (maxi != -1 && (/*foooonf ||*/ anchors[maxi].second)) {
			// if (maxi == 313771) {
			// 	foooonf = 1;
			// 	eprn("woohoo {}", paths.size());
			// }
			paths.back().push_front(maxi);
			anchors[maxi].second = false;
			maxi = prev[maxi];
		}

		int qlo = anchors[paths.back().front()].first.query_start, 
			qhi = anchors[paths.back().back()].first.query_end;
		int rlo = anchors[paths.back().front()].first.ref_start,   
			rhi = anchors[paths.back().back()].first.ref_end;

// ||> 470,967..478,485 -> 507,252..514,823;   7,518   7,571; e=10.2 g= 2.0 m= 8.3
		if (min(rhi - rlo, qhi - qlo) >= (1 - MAX_ERROR) * MIN_READ_SIZE) {
			dprn("-- chain: (len:{}) {}..{} --> {}..{} [dp={} exp_size={}]", abs(qhi-qlo), 
			qlo, qhi, rlo, rhi, dp[m.second], size_so_far[m.second]);
		}

		// if (foooonf) {
		// 	eprn("path {}", paths.back().size());
		// 	for (auto i: paths.back()) {
		// 		dprn("=== [{:4}] {:6}..{:6} -> {:6}..{:6} -> {}", i, 
		// 			anchors[i].first.query_start, anchors[i].first.query_end,
		// 			anchors[i].first.ref_start, anchors[i].first.ref_end,
		// 			dp[i]
		// 		);
		// 	}
		// }
	}

	return paths;
}

// ||> 470,967..478,485 -> 507,252..514,823;   7,518   7,571; e=10.2 g= 2.0 m= 8.3
/******************************************************************************/

vector<Hit> fast_align(const string &query, const string &ref, int kmer_size)
{
	auto T = cur_time();
	dprn("-- aligning query {:n} --> ref {:n}", query.size(), ref.size());

	// 1. Generate the list of hits (small anchors) inside the dot graph	
	auto anchors = generate_anchors(query, ref, kmer_size);

	// 2. Run DP on the anchors and collect all different anchors
	vector<Anchor> chains;
	auto chains_init = chain_anchors(anchors);	
	for (auto &chain: chains_init) {
		int qlo = anchors[chain.front()].first.query_start, 
			qhi = anchors[chain.back()].first.query_end;
		int rlo = anchors[chain.front()].first.ref_start,   
			rhi = anchors[chain.back()].first.ref_end;

		if (min(rhi - rlo, qhi - qlo) < (1 - 0.5) * MIN_READ_SIZE)
			continue;

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
	}
	dprn(":: elapsed/dp = {}s", elapsed(T)); T=cur_time();
	// exit(0);

	// 3. Perform the full alignment
	dprn("-- initial chains {}", chains.size());
	dprn(":: elapsed/long = {}s", elapsed(T)); T=cur_time();
	vector<Hit> hits;

	auto query_ptr = make_shared<Sequence>("QRY", query);
	auto ref_ptr = make_shared<Sequence>("REF", ref);
	for (auto &ch: chains) {
		auto hit = Hit {
			query_ptr, ch.query_start, ch.query_end,
			ref_ptr, ch.ref_start, ch.ref_end,
			0, "", "", {}
		};
		//if (false) {
		hit.aln = Alignment::from_anchors(query, ref, ch.query_kmers, ch.ref_kmers, 0);
		// hit.comment += fmt::format("ORIG:{}..{}->{}..{}",
		// 	ch.query_start, ch.query_end,
		// 	ch.ref_start, ch.ref_end);
		hit.query_start = hit.aln.start_a;
		hit.ref_start = hit.aln.start_b;
		hit.query_end = hit.aln.end_a;
		hit.ref_end = hit.aln.end_b;
		hits.push_back(hit);
		//}
	}
	dprn(":: elapsed/alignment = {}s", elapsed(T)); T=cur_time();

	return hits;
}


/******************************************************************************/

 // 9494..8508 [+] 21 -> 19385 atatcacagtgggtgtacacc
 // 9,515 --> ref 8,529

// 478214..514552 [+] 21 -> 438446 -> atatcacagtgggtgtacacc ~ 313771



void test(int, char** argv)
{
	FastaReference fr("data/hg19/hg19.fa");
	string s;

	Hit orig = Hit::from_bed(

		// "chr11	3613666	3623181	chr3	125420867	125429396	align_both/0001/both008387	0	+	-"
		// "chr10	98319847	98321003	chr10	126555985	126557221	align_both/0001/both006029	0	+	+	0	0	0"
		"chr1	88000	121417	chr1	235525	267707	align_both/0012/both060569	0	0	+	+"
	);
	Hit h = Hit::from_bed(
		// "chr11	3144946	3696247	chr3	125374462	125935440	0	0	+	-	560978	0	OK;;"
		// "chr10	98308557	98331750	chr10	126545044	126567993	0	0	+	+	23193	0			OK"
		"chr1	50861	154752	chr1	211516	283810	0	0	+	+	103891"
	);
	eprn("link: http://humanparalogy.gs.washington.edu/build37/{}", orig.name);
	eprn("RC: {} {}", h.query->is_rc, h.ref->is_rc);
	// h = orig;

	auto q = fr.get_sequence(h.query->name, h.query_start, h.query_end);
	auto r = fr.get_sequence(h.ref->name, h.ref_start, h.ref_end);
	if (h.ref->is_rc) r = rc(r);
	XX=q;
	
	const int k = 11;

	eprn("{} {}", string(60, '*'), k);
	eprn("Q -> {}...{}\nR -> {}...{}", q.substr(0, 20), q.substr(q.size() - 20),
		r.substr(0, 20), r.substr(r.size() - 20));
	eprn("{} {}", string(60, '*'), k);

	auto T = cur_time();
	auto hits = fast_align(q, r, k);
	sort(hits.begin(), hits.end(), [](Hit a, Hit b){ return a.query_start < b.query_start; });
	eprn("{} {}", string(60, '*'), k);
	for (auto &hit: hits) {
		eprn("||> {:7n}..{:7n} -> {:7n}..{:7n}; {:7n} {:7n}; e={:4.1f} g={:4.1f} m={:4.1f}", 
			hit.query_start, hit.query_end,
			hit.ref_start, hit.ref_end,
			hit.query_end - hit.query_start, hit.ref_end - hit.ref_start,
			hit.aln.error.error(), hit.aln.error.gap_error(), hit.aln.error.mis_error(),
			hit.comment
		);
		// eprn("{}", hit.aln.print(80));
	}

	int oqs = orig.query_start - h.query_start, 
		 oqe = orig.query_end - h.query_start,
		 ors = !orig.ref->is_rc ? orig.ref_start - h.ref_start : (h.ref_end - h.ref_start) - (orig.ref_end   - h.ref_start), 
		 ore = !orig.ref->is_rc ? orig.ref_end   - h.ref_start : (h.ref_end - h.ref_start) - (orig.ref_start - h.ref_start);

	eprn("{} {}", string(60, '*'), k);
	eprn("Looking for {:n}..{:n} -> {:n}..{:n}", oqs, oqe, ors, ore);
	eprn("Q -> {}...{}\nR -> {}...{}", q.substr(oqs, 20), q.substr(oqe - 20, 20),
		r.substr(ors, 20), r.substr(ore - 20, 20));
	
	//eprn("total {}s\n", elapsed(T));

	// string q = "ATCCTTGAAGCGCCCCCAAGGGCATCTTCTCAAAGTTGGATGTGTGCATTTTCCTGAGAGGAAAGCTTTCCCACATTATACAGCTTCTGAAAGGGTTGCTTGACCCACAGATGTGAAGCTGAGGCTGAAGGAGACTGATGTGGTTTCTCCTCAGTTTCTCTGTGTGGCACCAGGTGGCAGCAGAGGTCAGCAAGGCAAACCCGAGCCCAGGGATGCGGGGTGGGGGCAGGTACATCCTCTCTTGAGCTACAGCAGATTAACTCTGTTCTGTTTCATTGTGGTTGTTTAGTTTGCGTTTTTTTTTCTCCAACTTTGTGCTTCATCGGGAAAAGCTTTGGATCACAATTCCCAGTGCTGAAGAAAAGGCCAAACTCTGGAAAAAATTTGAATATTTTGAGCCAAATGTGAGGACCACAACCTGTGAGAACGGAAAATAAATCCTGGGACCCCAGACTCACTAAGCCAAAGGGAAAAGCCAAGCTGGGAACTGGCTTATGCAAACCTGCTTCCCATCTGGTTCCTAAATAAGATAGCTATTACACAAAGACAAAAAAGCTACATCCCTGCCTCTACCTCCATCGCATGCAAAATGTGTATTCAGTGAACGCTGACCAAAGACAGAAGAATGCAACCATTTGCCTCTGATTTACCCACACCCATTTTTTCCACTTCTTCCCCTTTCCCCAATACCCGCACTTTTCCCCTTTACTTACTGAGGTCCCCAGACAACCTTTGGGAAAAGCACGGACCACAGTTTTTCCTGTGGTTCTCTGTTCTTTTCTCAGGTGTGTCCTTAACCTTGCAAATAGATTTCTTGAAATGATTGAGACTCACCTTGGTTGTGTTCTTTGATTAGTGCCTGTGACGCAGCTTCAGGAGGTCCTGAGAACGTGTGCACAGTTTAGTCGGCAGAAACTTAGGGAAATGTAAGACCACCATCAGCACATAGGAGTTCTGCATTGGTTTGGTCTGCATTGGTTTGGTCTGGAAGGAGGAAAATTCAAAGTAATGGGGCTTACAGGTCATAGATAGATTCAAAGATTTTCTGATTGTCAATTGGTTGAAAGAATTATTATCTACAGACCTGCTATCAATAGAAAGGAGAGTCTGGGTTAAGATAAGAGACTGTGGAGACC";
	// string r = "ATCCTTGAAGCGCCCCCAAGGGCATCTTCTCAAAGTTGGATGTGTGCATTTTCCTGAGAGGAAAGCTTTCCCACATTATTCAGCTTCTGAAAGGGTTGCTTGACCCACAGATGTGAAGCTGAGGCTGAAGGAGACTGATGTGGTTTCTCCTCAGTTTCTCTGTGCGGCACCAGGTGGCAGCAGAGGTCAGCAAGGCAAACCCGAGCCCGGGGATGCGGGGTGGGGGCAGCTACGTCCTCTCTTGAGCTACAGCAGATTCACTCTGTTCTGTTTCATTGTTGCTTAGTTTGCGTTTTGTTTCTCCAACTTTGTGCCTCATCAGGAAAAGCTTTGGATCACAATTCCCAGTGCTGAAGAAAAGGCCAAACTCTGGAAAAAATTTTGAATATTTTGAGCCAAATGTGAGGACCACAACCTGTGAGAACGGAAAATAAATCCTGGGACCCCAGACTCACTAAGCCAAAGGGAAAAGCCAAGCTGGGAACTGGCTTATGCAAACCTGCTTCCCATCTGGTTCCTAAATAAGATAGCTATTACACAAAGATAAAAAAGCTACATCCCTGCCTCTACCTCCCTCGCATGTAAAATGTGTATTCAGTGAACACTGACCAAAGACAGAAGAATGCAACCATTTGCCTCTGATTTACCCACACCCATTTTTTCCACTTCTTCCCCTTTCCCCAATACCCGCACTTTTCCCCTTTACTTACTGAGGCCCCCAGACAATCTTTGGGAAAAGCACGGACCACAGTTTTTCCTGTGGTTCTCTGTTCTTTTCTCAGGTGTGTCCTTAACCTTGCAAATAGATTTCTTGAAATGATTGACACTCACCTTGGTTGTGTTCTTTGATCAGCGCCTGTGACGCAGCTTCAGGAGGTCCTGAGAACGTGTGCACAGTTTAGTCGGCAGAAACTTAGGGAAACGTAAGACCACCATCAGTACGTAGGAGTTGTGCATTGGTTTGGTCTGGAAGGAGGAAAATTCAAAGTAATGGGGCTTACAGGTCATAGATAGATTCAAAGATTTTCTGATTGTCAATTGATTGAAAGAATTATTATCTACAGACCTGCTATCAATAGAAAGGAGAGTCTGAGTTAAGATAAGAGACTGTGGAGACC";

	// int qs, qe, rs, re;
	// string q = fr.get_sequence("chr21", qs=10596500-1000, qe=10598323+1000);
	// string r = fr.get_sequence("chr19", rs=37763913-1000, re=37765775+1000);

	// // 1836 1862

	// // string q = fr.get_sequence("chr21", qs=10596477-1000, qe=10598323+1000);
	// // string r = fr.get_sequence("chr19", rs=37761372-1000, re=37763248+1000);
	// dprn("{}..{} --> {}..{}", qs, qe, rs, re);
	
	// shared_ptr<Index> query_hash = make_shared<Index>(make_shared<Sequence>("qry", q), 12);
	// shared_ptr<Index> ref_hash = make_shared<Index>(make_shared<Sequence>("ref", r), 12);

	// Tree tree;
	// auto hi = search(0, query_hash, ref_hash, tree, false, 1800, false);
	// 	for (auto &pp: hi) {
	// 	dprn("{}..{} -> {}..{} // {} ~ {}",
	// 		pp.query_start, pp.query_end,
	// 		pp.ref_start, pp.ref_end,
	// 		abs(pp.query_start-pp.query_end),
	// 		abs(pp.ref_start-pp.ref_end)
	// 	);
	// }
	// dprn("************************************************************");

	// Tree tree2;
	// for (int qi = 0; qi < query_hash->minimizers.size(); qi++) {
	// 	auto &qm = query_hash->minimizers[qi];
	// 	if (qm.hash.status != Hash::Status::HAS_UPPERCASE) 
	// 		continue; 
	// 	auto hi = search(qi, query_hash, ref_hash, tree2, false);
	// 	for (auto &pp: hi) {
	// 		dprn("{}..{} -> {}..{} // {} ~ {}",
	// 			pp.query_start, pp.query_end,
	// 			pp.ref_start, pp.ref_end,
	// 			abs(pp.query_start-pp.query_end),
	// 			abs(pp.ref_start-pp.ref_end)
	// 		);
	// 	}
}
