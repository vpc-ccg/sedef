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
	vector<Hit> anchors;

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
					if (query[q + len] == 'N' || ref[r + len] == 'N')
						break;
					has_u |= (isupper(query[q + len]) || isupper(ref[r + len]));
					if (toupper(query[q + len]) != toupper(ref[r + len]))
						break;
				}
				if (len >= kmer_size / 2) {
					anchors.push_back(Hit{
						nullptr, q, q + len, 
						nullptr, r, r + len, 
						has_u
					});
					slide[d] = q + len;
				}
			} else {
				assert(slide[d] >= q + kmer_size);
				assert(d - off + slide[d] >= r + kmer_size);
			}
		}
	}
	
	for (int i = 1; i < anchors.size(); i++) {
		assert(tie(anchors[i - 1].query_start, anchors[i - 1].ref_start) <= 
			tie(anchors[i].query_start, anchors[i].ref_start));
	}
	return anchors; 
}

/******************************************************************************/

/*TODO*/ string XX,YY;

auto chain_anchors(vector<Hit> &anchors)
{
	const double ratio = 4; // ratio: match to error
	const int max_gap = MAX_ERROR * MIN_READ_SIZE;

	auto T = cur_time();

	struct Coor {
		pair<int, int> x; 
		int score, pos;
		bool operator<(const Coor &a) const { return x < a.x; }
	};
	vector<Coor> xs, ys; // QUERY; pos in ANCHORS
	xs.reserve(2 * anchors.size());
	ys.reserve(anchors.size());
	
	int max_q = 0, max_r = 0, l = 0;
	for (int i = 0; i < anchors.size(); i++) {
		l++;
		auto &a = anchors[i];
		xs.push_back({{a.query_start, i}, numeric_limits<int>::min(), i});
		xs.push_back({{a.query_end, i}, numeric_limits<int>::min(), i});
		ys.push_back({{a.ref_end - 1, i}, numeric_limits<int>::min(), i}); 
		
		assert(a.query_start < a.query_end);
		max_q = max(max_q, a.query_end);
		max_r = max(max_r, a.ref_end);
	}
	// for (auto a: anchors) 
	// 	dprn("!== {:6}..{:6} -> {:6}..{:6}", a.query_start, a.query_end, a.ref_start, a.ref_end);

	dprn("-- anchors to dp: {:n}", l);

	sort(xs.begin(), xs.end());
	SegmentTree<Coor> tree(ys); 

	vector<int> prev(anchors.size(), -1);
	vector<pair<int, int>> dp(anchors.size());
	for (int i = 0; i < dp.size(); i++)
		dp[i] = {0, i};
	// vector<int> size_so_far(anchors.size(), 0);

	// set<pair<int, int>, greater<pair<int, int>>> maxes;
	int deactivate_bound = 0;

	// dprn(">>>> init {}", elapsed(T)); T=cur_time();
	for (auto &x: xs) {
		int i = x.x.second;
		auto &a = anchors[i];
		if (x.x.first == a.query_start) {
			while (deactivate_bound < (&x - &xs[0])) {
				int t = xs[deactivate_bound].x.second; // index
				if (xs[deactivate_bound].x.first == anchors[t].query_end) { // end point
					if (a.query_start - anchors[t].query_end <= max_gap)
						break;
					tree.deactivate({anchors[t].ref_end - 1, t});
				}
				deactivate_bound++;
			}

			int w = ratio * (a.query_end - a.query_start);
			int j = tree.rmq({a.ref_start - max_gap, 0}, 
				             {a.ref_start - 1, anchors.size()});
			if (j != -1 && ys[j].score != numeric_limits<int>::min()) {
				j = ys[j].pos;
				auto &p = anchors[j];
				assert(a.query_start >= p.query_end);
				assert(a.ref_start >= p.ref_end);
				int gap = (a.query_start - p.query_end + a.ref_start - p.ref_end);
				if (w + dp[j].first - gap > 0) {
					dp[i].first = w + dp[j].first - gap;
					prev[i] = j;
					// size_so_far[i] = size_so_far[j] + (a.query_end - a.query_start) + 
					// 	max(a.query_start - p.query_end, a.ref_start - p.ref_end);
				} else {
					dp[i].first = w;
					// size_so_far[i] = (a.query_end - a.query_start);
				}
			} else {
				dp[i].first = w;
				// size_so_far[i] = (a.query_end - a.query_start);
			}
			//if (dp[i] >= (MIN_READ_SIZE / 2) * (ratio * (1 - MAX_ERROR) - MAX_ERROR))
			// maxes.insert({dp[i].first, i});
		} else {
			int gap = (max_q + 1 - a.query_end + max_r + 1 - a.ref_end);
			tree.activate({a.ref_end - 1, i}, dp[i].first - gap);
		}
	}
	sort(dp.begin(), dp.end(), greater<pair<int, int>>());
	// dprn(">>>> search {}", elapsed(T)); T=cur_time();

	vector<int> path; path.reserve(anchors.size());
	vector<pair<int, bool>> boundaries {{0, 0}};
	vector<char> used(anchors.size(), 0);
	for (auto &m: dp) {
		int maxi = m.second;
		if (used[maxi])
			continue;
		bool has_u = 0;
		while (maxi != -1 && !used[maxi]) {
			path.push_back(maxi);
			// paths.back().first.push_front(maxi);
			has_u |= anchors[maxi].jaccard;
			used[maxi] = true;
			maxi = prev[maxi];
		}
		for (int ai = path.size() - 2; ai >= boundaries.back().first; ai--) {
			assert(anchors[path[ai + 1]] < anchors[path[ai]]);
		}
		boundaries.push_back({path.size(), has_u});
	}
	// dprn(">>>> recon {}", elapsed(T)); T=cur_time();
	return make_pair(path, boundaries);
}

/******************************************************************************/

vector<Hit> fast_align(const string &query, const string &ref, int kmer_size)
{
	auto T = cur_time();
	dprn("-- aligning query {:n} --> ref {:n}", query.size(), ref.size());
	auto query_ptr = make_shared<Sequence>("QRY", query);
	auto ref_ptr = make_shared<Sequence>("REF", ref);

	/// 1. Generate the list of hits (small anchors) inside the dot graph	
	auto anchors = generate_anchors(query, ref, kmer_size);

	/// 2. Run DP on the anchors and collect all different anchors
	vector<Hit> hits;
	auto chains_init = chain_anchors(anchors);
	auto &bounds = chains_init.second;	
	auto &chain = chains_init.first;
	for (int bi = 1; bi < bounds.size(); bi++) {
		bool has_u = bounds[bi].second;
		int be = bounds[bi].first;
		int bs = bounds[bi - 1].first;

		int qlo = anchors[chain[be - 1]].query_start, 
			qhi = anchors[chain[bs]].query_end;
		int rlo = anchors[chain[be - 1]].ref_start,   
			rhi = anchors[chain[bs]].ref_end;

		// check error
		int span = max(rhi - rlo, qhi - qlo);
		if (!(has_u && span >= 90) && !(span >= MIN_READ_SIZE * (1 - MAX_ERROR)))
			continue;

		assert(qhi <= query.size());
		assert(rhi <= ref.size());

		Hit a { query_ptr, qlo, qhi, ref_ptr, rlo, rhi, 0, "", "", {}, {} };
		for (int bi = be - 1; bi >= bs; bi--) {
			int ai = chain[bi];
			a.guides.first.push_back({anchors[ai].query_start, anchors[ai].query_end});
			a.guides.second.push_back({anchors[ai].ref_start, anchors[ai].ref_end});
		}
		hits.push_back(a);
	}
	dprn(":: elapsed/dp = {}s", elapsed(T)); T=cur_time();

	/// 3. Perform the full alignment
	for (auto &hit: hits) {
		hit.aln = Alignment::from_anchors(query, ref, hit.guides.first, hit.guides.second, 0);
		hit.query_start = hit.aln.start_a;
		hit.ref_start = hit.aln.start_b;
		hit.query_end = hit.aln.end_a;
		hit.ref_end = hit.aln.end_b;
	}
	dprn(":: elapsed/alignment = {}s", elapsed(T)); T=cur_time();

	/// 3. Refine these chains
	void refine_chains(vector<Hit> &anchors);
	refine_chains(hits);
	dprn(":: elapsed/refinement = {}s", elapsed(T)); T=cur_time();

	return hits;
}


/******************************************************************************/

 // 9494..8508 [+] 21 -> 19385 atatcacagtgggtgtacacc
 // 9,515 --> ref 8,529

// 478214..514552 [+] 21 -> 438446 -> atatcacagtgggtgtacacc ~ 313771

#include <fstream>
void test(int, char** argv)
{
	// auto x = align(
	// 	"CAAGAGAATTAAATGGGTTATTGATTAAAAA",
	// 	"TTTTTTCAAGAGAATTAAATCATTTCTTGATTA"
	// );
	// eprn("{}", x.print());
	// x.trim_back();
	// eprn("{}", x.print());
	// x.trim_front();
	// eprn("{}", x.print());
	// exit(0);

	FastaReference fr("data/hg19/hg19.fa");
	string s, s2, sl;
	const int k = 11;

	ifstream MISS("out/misses.txt");
	int line = 0;
	while (getline(MISS, sl)) {
		auto TT = cur_time();
		auto ssl = split(sl, ' ');
		if (ssl[1] != "align_both/0013/both067680") continue;

		ifstream fin("out/bucket_0000_temp_diff.bed");
		s2 = "";
		while (getline(fin, s)) {
			auto ss = split(s, '\t');
			if (ss[6] == ssl[1]) {
				for (int x = 10; x < ss.size(); x++) s2 += ss[x] + "\t";
				break;
			}
		}
		
		eprn("Hit http://humanparalogy.gs.washington.edu/build37/{} @ {}", ssl[1], ssl[0]);
		if (s2 == "") {
			eprn("-- complete miss");
			continue;
		}
		Hit orig = Hit::from_bed(
			// "chr11	3613666	3623181	chr3	125420867	125429396	align_both/0001/both008387	0	+	-"
			// "chr10	98319847	98321003	chr10	126555985	126557221	align_both/0001/both006029	0	+	+	0	0	0"
			// "chr1	88000	121417	chr1	235525	267707	align_both/0012/both060569	0	0	+	+"
			s
		);
		Hit h = Hit::from_bed(
			// "chr11	3144946	3696247	chr3	125374462	125935440	0	0	+	-	560978	0	OK;;"
			// "chr10	98308557	98331750	chr10	126545044	126567993	0	0	+	+	23193	0			OK"
			// "chr1	50861	154752	chr1	211516	283810	0	0	+	+	103891"
			s2
		);
		eprn("{}{} {}...", "+-"[orig.query->is_rc], "+-"[orig.ref->is_rc], orig.to_bed(false).substr(0, 75));
		eprn("{}{} {}...", "+-"[h.query->is_rc], "+-"[h.ref->is_rc], h.to_bed(false).substr(0, 50));

		auto q = fr.get_sequence(h.query->name, h.query_start, h.query_end);
		auto r = fr.get_sequence(h.ref->name, h.ref_start, h.ref_end);
		if (h.ref->is_rc) r = rc(r);

		eprn("{} {}", string(60, '*'), k);
		int oqs = orig.query_start - h.query_start, 
			oqe = orig.query_end - h.query_start,
			ors = !orig.ref->is_rc ? orig.ref_start - h.ref_start : (h.ref_end - h.ref_start) - (orig.ref_end   - h.ref_start), 
			ore = !orig.ref->is_rc ? orig.ref_end   - h.ref_start : (h.ref_end - h.ref_start) - (orig.ref_start - h.ref_start);

		vector<int> oqcov(oqe - oqs, 0), 
					orcov(ore - ors, 0);
		auto T = cur_time();
		auto hits = fast_align(q, r, k);
		sort(hits.begin(), hits.end(), [](Hit a, Hit b){ return a.query_start < b.query_start; });
		eprn("{} {}", string(60, '*'), k);
		for (auto &hit: hits) {
			// eprn("||> {:7n}..{:7n} -> {:7n}..{:7n}; {:7n} {:7n}; e={:4.1f} g={:4.1f} m={:4.1f}", 
			// 	hit.query_start, hit.query_end,
			// 	hit.ref_start, hit.ref_end,
			// 	hit.query_end - hit.query_start, hit.ref_end - hit.ref_start,
			// 	hit.aln.error.error(), hit.aln.error.gap_error(), hit.aln.error.mis_error(),
			// 	hit.comment
			// );

			int over_qs = max(hit.query_start, oqs),
				over_qe = min(hit.query_end, oqe);
			int over_rs = max(hit.ref_start, ors),
				over_re = min(hit.ref_end, ore);

			if (over_qs <= over_qe && over_rs <= over_re) {
				for (int i = over_qs; i < over_qe; i++)
					oqcov[i - oqs]++;
				for (int i = over_rs; i < over_re; i++)
					orcov[i - ors]++;
			}

			// eprn("{}", hit.aln.print(80));
		}

		eprn("{} {}", string(60, '*'), k);

		eprn("-- looking for {:n}..{:n} -> {:n}..{:n}", oqs, oqe, ors, ore);
		eprn("   Q -> {}...{}\n   R -> {}...{}", q.substr(oqs, 20), q.substr(oqe - 20, 20),
			r.substr(ors, 20), r.substr(ore - 20, 20));

		int p=0;
		for (int i = 0; i < oqcov.size(); ) {
			int j;
			for (j = i + 1; j < oqcov.size() && oqcov[i] == oqcov[j]; j++);
			eprnn("   |> {:6}..{:6}~>{:2.0f}%/{}={:<2} ", i + oqs, j + oqs, 
				pct(j-i,oqcov.size()), j-i, oqcov[i]);
			if (oqcov[i]) p+=pct(j-i,oqcov.size());
			i = j;
		} eprn("=== {}", p); int p1=p; p = 0;
		for (int i = 0; i < orcov.size(); ) {
			int j;
			for (j = i + 1; j < orcov.size() && orcov[i] == orcov[j]; j++);
			eprnn("   |> {:6}..{:6}~>{:2.0f}%/{}={:<2} ", i + ors, j + ors, pct(j-i,orcov.size()), 
				j-i, orcov[i]);
			if (orcov[i]) p+=pct(j-i,orcov.size());
			i = j;
		} eprn("=== {}", p);

		eprn("line {:3} {} @ {} : {:5} hits, {:3} {:3} --> {:3.1f} s\n", line++, ssl[1], ssl[0],
			hits.size(), p1, p, elapsed(TT));
		//if (p1 < 95 || p < 95)
			// cin.get();
	}
	//eprn("total {}s\n", elapsed(T));

	// -- looking for 20,321..22,929 -> 19,569..23,541

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
