/// 786

#include <fstream>

#include "common.h"
#include "fasta.h"
#include "align.h"
#include "search.h"
#include "segment.h"
#include "chain.h"
#include "refine.h"

using namespace std; 

/******************************************************************************/

const int MIN_UPPERCASE_MATCH = 90;
const int MAX_CHAIN_GAP = MAX_ERROR * MIN_READ_SIZE;

const int MATCH_CHAIN_SCORE = 4;

/******************************************************************************/

vector<Anchor> generate_anchors(const string &query, const string &ref, const int kmer_size)
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
	vector<Anchor> anchors;
	dprn("got {} hashes", ref_hashes.size());

	// vector<int> qq;
	// for(auto &q:ref_hashes) {qq.push_back(q.second.size());}
	// 	sort(qq.begin(), qq.end());
	// eprn("{}", qq.back());qq.pop_back();
	// eprn("{}", qq.back());qq.pop_back();
	// eprn("{}", qq.back());

	last_n = -kmer_size, h = 0;
	int w = 0;
	for (int i = 0; i < query.size(); i++) {
		if (query[i] == 'N') 
			last_n = i;
		h = ((h << 2) | hash_dna(query[i])) & MASK; 
		if (i < kmer_size - 1) 
			continue;
		if (last_n >= (i - kmer_size + 1)) 
			continue;
		
		auto it = ref_hashes.find(h);
		if (it == ref_hashes.end()) // || it->second.size() >= 1000)
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
					if (query[q + len] == 'N' || ref[r + len] == 'N') {
						assert(len >= kmer_size);
						break;
					}
					if (toupper(query[q + len]) != toupper(ref[r + len]))
						break;
					has_u |= (isupper(query[q + len]) || isupper(ref[r + len]));
				}
				if (len >= kmer_size) {
					if (anchors.size() >= (1<<20) && anchors.size() == anchors.capacity()) {
						anchors.reserve(anchors.size() * 1.5);
					}
					anchors.emplace_back(Anchor{q, r, len, has_u});
					slide[d] = q + len;
				}
			} else {
				assert(slide[d] >= q + kmer_size);
				assert(d - off + slide[d] >= r + kmer_size);
			}
		}
	}
	dprn("{:n} {:n}", anchors.size(), anchors.capacity());
	// for (int i = 1; i < anchors.size(); i++) {
	// 	assert(tie(anchors[i - 1].query_start, anchors[i - 1].ref_start) <= 
	// 		tie(anchors[i].query_start, anchors[i].ref_start));
	// }
	return anchors; 
}

auto chain_anchors(vector<Anchor> &anchors)
{
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
		xs.push_back({{a.q, i}, SegmentTree<Coor>::MIN, i});
		xs.push_back({{a.q + a.l, i}, SegmentTree<Coor>::MIN, i});
		ys.push_back({{a.r + a.l - 1, i}, SegmentTree<Coor>::MIN, i}); 
		
		assert(a.q < a.q + a.l);
		max_q = max(max_q, a.q + a.l);
		max_r = max(max_r, a.r + a.l);
	}
	// for (auto a: anchors) 
	// 	dprn("!== {:6}..{:6} -> {:6}..{:6}", a.query_start, a.query_end, a.ref_start, a.ref_end);
	dprn("-- anchors to dp: {:n}", l);

	sort(xs.begin(), xs.end());
	SegmentTree<Coor> tree(ys); 

	vector<int> prev(anchors.size(), -1);
	vector<pair<int, int>> dp(anchors.size());
	for (int i = 0; i < dp.size(); i++) {
		dp[i] = {0, i};
	}
	int deactivate_bound = 0;
	for (auto &x: xs) {
		int i = x.x.second;
		auto &a = anchors[i];
		if (x.x.first == a.q) {
			while (deactivate_bound < (&x - &xs[0])) {
				int t = xs[deactivate_bound].x.second; // index
				if (xs[deactivate_bound].x.first == anchors[t].q + anchors[t].l) { // end point
					if (a.q - anchors[t].q + anchors[t].l <= MAX_CHAIN_GAP)
						break;
					tree.deactivate({anchors[t].r + anchors[t].l - 1, t});
				}
				deactivate_bound++;
			}

			int w = MATCH_CHAIN_SCORE * (a.q + a.l - a.q);
			int j = tree.rmq({a.r - MAX_CHAIN_GAP, 0}, 
				             {a.r - 1, anchors.size()});
			if (j != -1 && ys[j].score != SegmentTree<Coor>::MIN) {
				j = ys[j].pos;
				auto &p = anchors[j];
				assert(a.q >= p.q + p.l);
				assert(a.r >= p.r + p.l);
				int gap = (a.q - p.q + p.l + a.r - p.r + p.l);
				if (w + dp[j].first - gap > 0) {
					dp[i].first = w + dp[j].first - gap;
					prev[i] = j;
				} else {
					dp[i].first = w;
				}
			} else {
				dp[i].first = w;
			}
		} else {
			int gap = (max_q + 1 - a.q + a.l + max_r + 1 - a.r + a.l);
			tree.activate({a.r + a.l - 1, i}, dp[i].first - gap);
		}
	}
	sort(dp.begin(), dp.end(), greater<pair<int, int>>());

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
			has_u |= anchors[maxi].has_u;
			used[maxi] = true;
			maxi = prev[maxi];
		}
		// for (int ai = path.size() - 2; ai >= boundaries.back().first; ai--) {
		// 	assert(anchors[path[ai + 1]] < anchors[path[ai]]);
		// }
		boundaries.push_back({path.size(), has_u});
	}
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
	dprn("-- got {} anchors in {} s", anchors.size(), elapsed(T)); T=cur_time();

	/// 2. Run DP on the anchors and collect all different anchors
	vector<Hit> hits;
	vector<vector<int>> guides;
	auto chains_init = chain_anchors(anchors);
	auto &bounds = chains_init.second;	
	auto &chain = chains_init.first;
	for (int bi = 1; bi < bounds.size(); bi++) {
		bool has_u = bounds[bi].second;
		int be = bounds[bi].first;
		int bs = bounds[bi - 1].first;

		int qlo = anchors[chain[be - 1]].q,
			qhi = anchors[chain[bs]].q + anchors[chain[bs]].l;
		int rlo = anchors[chain[be - 1]].r,   
			rhi = anchors[chain[bs]].r + anchors[chain[bs]].l;

		// check error
		int span = max(rhi - rlo, qhi - qlo);
		if (!(has_u && span >= MIN_UPPERCASE_MATCH) && !(span >= MIN_READ_SIZE * (1 - MAX_ERROR)))
			continue;

		assert(qhi <= query.size());
		assert(rhi <= ref.size());

		Hit a { query_ptr, qlo, qhi, ref_ptr, rlo, rhi };
		guides.push_back(vector<int>());
		for (int bi = be - 1; bi >= bs; bi--) {
			guides.back().push_back(chain[bi]);
		}
		hits.push_back(a);
	}
	dprn(":: elapsed/dp = {}s", elapsed(T)); T=cur_time();
	
	/// 3. Perform the full alignment
	for (auto &hit: hits) {
		hit.aln = Alignment(query, ref, anchors, guides[&hit - &hits[0]]);
		update_from_alignment(hit);
	}
	dprn(":: elapsed/alignment = {}s", elapsed(T)); T=cur_time();

	/// 3. Refine these chains
	refine_chains(hits, query, ref);
	dprn(":: elapsed/refinement = {}s", elapsed(T)); T=cur_time();

	return hits;
}


/******************************************************************************/

void test2()
{
// 1556 out of 1559 (99.8, len 115473..123766)      chr16	32019249	32203896	chr16	33746594	33929593			+	+	184647	0			OK;;;
// 1557 out of 1559 (99.9, len 184647..182999)      chr16	5128906	5451300	chr4	3946932	4274093			+	-327161	0			OK;;;

	const int k=11;
	FastaReference fr("data/hg19/hg19.fa");

	auto TT = cur_time();
	Hit h = Hit::from_bed(
		"chrX	35014	3024029	chrY	0	2984783			+	+	2989015	0	;;"
	);
	eprn("{}{} {}...", "+-"[h.query->is_rc], "+-"[h.ref->is_rc], h.to_bed(false).substr(0, 50));

	auto q = fr.get_sequence(h.query->name, h.query_start, &h.query_end);
	auto r = fr.get_sequence(h.ref->name, h.ref_start, &h.ref_end);
	if (h.ref->is_rc) r = rc(r);
	assert(r.size() == h.ref_end - h.ref_start);
	assert(q.size() == h.query_end - h.query_start);

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
			hit.aln.total_error(), hit.aln.gap_error(), hit.aln.mismatch_error(),
			hit.comment
		);
	}
	eprn("done in {} s", elapsed(TT));
}

void test(int, char** argv)
{
	// test2();
	// exit(0);
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
		auto ssl = vector<string>{"miss", sl}; // split(sl, ' ');
		// if (ssl[1] != "align_both/0019/both097366") continue;

		ifstream fin("out/ww.bed");
		bool ok=0;
		while (getline(fin, s)) {
			auto ss = split(s, '\t');
			if (ss[6] == ssl[1]) {
				s2 = "";
				for (int x = 10; x < ss.size(); x++) s2 += ss[x] + "\t";
				if (ss[0] == ss[10] && ss[3] == ss[13]) {
					ok=1;
					break;
				}
			}
		}
		if (!ok) continue;
		
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

		if (h.query_start >= orig.query_start) {
			int add = h.query_start - orig.query_start;
			h.query_start -= add;
			eprn(">>>>>> PADDED Q_L {}", add);
		}
		if (h.ref_start >= orig.ref_start) {
			int add = h.ref_start - orig.ref_start;
			h.ref_start -= add;
			eprn(">>>>>> PADDED R_L {}", add);
		}

		if (h.query_end < orig.query_end) {
			int add = orig.query_end - h.query_end;
			h.query_end += add;
			eprn(">>>>>> PADDED Q_R {}", add);
		}
		if (h.ref_end < orig.ref_end) {
			int add = orig.ref_end - h.ref_end;
			h.ref_end += add;
			eprn(">>>>>> PADDED R_R {}", add);
		}
		eprn("{}{} {}...", "+-"[orig.query->is_rc], "+-"[orig.ref->is_rc], orig.to_bed(false).substr(0, 75));
		eprn("{}{} {}...", "+-"[h.query->is_rc], "+-"[h.ref->is_rc], h.to_bed(false).substr(0, 50));

		auto ooq = fr.get_sequence(orig.query->name, orig.query_start, &orig.query_end);
		auto oor = fr.get_sequence(orig.ref->name, orig.ref_start, &orig.ref_end);
		if (orig.ref->is_rc) oor = rc(oor);

		// eprn("   Q -> {}...{}\n   R -> {}...{}", 
		// 	ooq.substr(0, 20), ooq.substr(ooq.size() - 20, 20),
		// 	oor.substr(0, 20), oor.substr(oor.size() - 20, 20));
		// eprn("{} {} {} {}", h.query->name, h.ref->name, h.ref->is_rc, orig.ref->is_rc);

		auto q = fr.get_sequence(h.query->name, h.query_start, &h.query_end);
		auto r = fr.get_sequence(h.ref->name, h.ref_start, &h.ref_end);
		if (h.ref->is_rc) r = rc(r);
		assert(r.size() == h.ref_end - h.ref_start);
		assert(q.size() == h.query_end - h.query_start);

		eprn("{} {}", string(60, '*'), k);
		int oqs = orig.query_start - h.query_start, 
			oqe = orig.query_end - h.query_start,
			ors = (!orig.ref->is_rc 
				? orig.ref_start - h.ref_start 
				: h.ref_end - orig.ref_end), 
			ore = (!orig.ref->is_rc 
				? orig.ref_end - h.ref_start 
				: h.ref_end - orig.ref_start);

// 511,319->520,991 q
// 391,794->542,894 r

		// eprn("{} {}->{} {}->{}", orig.ref->is_rc, orig.query_start, orig.query_end, orig.ref_start, orig.ref_end);
		// eprn("{} {}->{} {}->{}", h.ref->is_rc, h.query_start, h.query_end, h.ref_start, h.ref_end);

		int w = 0;
		eprn("--> {}", w=r.find(oor, w));
		assert(ooq == q.substr(oqs, oqe-oqs));
		assert(oor == r.substr(ors, ore-ors));

		assert(oqs>=0&&ors>=0);
		assert(oqe<=h.query_end-h.query_start);
		assert(ore<=h.ref_end-h.ref_start);

		vector<int> oqcov(oqe - oqs, 0), 
					orcov(ore - ors, 0);
		auto T = cur_time();
		auto hits = fast_align(q, r, k);
		sort(hits.begin(), hits.end(), [](Hit a, Hit b){ return a.query_start < b.query_start; });
		eprn("{} {}", string(60, '*'), k);
		for (auto &hit: hits) {
			eprn("||> {:7n}..{:7n} -> {:7n}..{:7n}; {:7n} {:7n}; e={:4.1f} g={:4.1f} m={:4.1f}", 
				hit.query_start, hit.query_end,
				hit.ref_start, hit.ref_end,
				hit.query_end - hit.query_start, hit.ref_end - hit.ref_start,
				hit.aln.total_error(), hit.aln.gap_error(), hit.aln.mismatch_error(),
				hit.comment
			);

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
		if (p1 < 80 || p < 80)
			cin.get();
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
