/// 786

#include "common.h"
#include "fasta.h"
#include "align.h"
#include "search.h"
#include "segment.h"
#include "chain.h"
#include "hit.h"

using namespace std; 

/******************************************************************************/

// 18,024..19,126 -> 17,524..19,268

// ||>  18,002.. 18,394 ->  17,505.. 17,892;     392     387; e= 9.9 g= 1.3 m= 8.7
// ||>  18,422.. 18,799 ->  18,232.. 18,612;     377     380; e= 9.5 g= 0.8 m= 8.7
// ||>  18,808.. 19,118 ->  18,950.. 19,260;     310     310; e= 5.2 g= 0.0 m= 5.2

void merge_alignments(Alignment &prev, Alignment &cur, const string &qstr, const string &rstr)
{
	assert(cur.start_a < prev.end_a || cur.start_b < prev.end_b);
	// eprn("merging... {}..{} to {}..{}", prev.start_a, prev.end_a, prev.start_b, prev.end_b);
	// eprn("         . {}..{} to {}..{}", cur.start_a, cur.end_a, cur.start_b, cur.end_b);

	assert(prev.end_a <= cur.end_a);
	assert(prev.end_b <= cur.end_b);
	
	int trim = prev.end_a - cur.start_a;
	// eprn("trim q: {}", trim);
	int q = 0, r = 0, i = 0;
	for (i = prev.alignment.size() - 1; i >= 0 && q < trim; i--) {
		if (prev.align_a[i] != '-') q++;
		if (prev.align_b[i] != '-') r++;
	}
	// eprn("i={}", i);
	prev.align_a = prev.align_a.substr(0, i+1);
	prev.alignment = prev.alignment.substr(0, i+1);
	prev.align_b = prev.align_b.substr(0, i+1);
	prev.end_a = prev.start_a + prev.a.size() - q;
	prev.end_b = prev.start_b + prev.b.size() - r;
	prev.a = prev.a.substr(0, prev.a.size() - q);
	prev.b = prev.b.substr(0, prev.b.size() - r);
	string X;
	// X="";for (auto c: prev.align_a) if (c!='-')X+=c; assert(X==prev.a); 
	// assert(prev.a==qstr.substr(prev.start_a, prev.end_a-prev.start_a));
	// X="";for (auto c: prev.align_b) if (c!='-')X+=c; assert(X==prev.b);
	// assert(prev.b==rstr.substr(prev.start_b, prev.end_b-prev.start_b));


	q = 0, r = 0, i = 0;
	for (i = 0; i < cur.alignment.size() && q < trim; i++) {
		if (cur.align_a[i] != '-') q++;
		if (cur.align_b[i] != '-') r++;
	}
	cur.align_a = cur.align_a.substr(i);
	cur.alignment = cur.alignment.substr(i);
	cur.align_b = cur.align_b.substr(i);
	cur.start_a += q;
	cur.start_b += r;
	cur.a = cur.a.substr(q);
	cur.b = cur.b.substr(r);
	// X="";for (auto c: cur.align_a) if (c!='-')X+=c; assert(X==cur.a);
		// assert(cur.a==qstr.substr(cur.start_a, cur.end_a-cur.start_a));
	// X="";for (auto c: cur.align_b) if (c!='-')X+=c; assert(X==cur.b);
		// assert(cur.b==rstr.substr(cur.start_b, cur.end_b-cur.start_b));


	trim = prev.end_b - cur.start_b;
	// eprn("trim r: {}", trim);
	q = 0, r = 0, i = 0;
	for (i = prev.alignment.size() - 1; i >= 0 && r < trim; i--) {
		if (prev.align_a[i] != '-') q++;
		if (prev.align_b[i] != '-') r++;
	}
	prev.align_a = prev.align_a.substr(0, i+1);
	prev.alignment = prev.alignment.substr(0, i+1);
	prev.align_b = prev.align_b.substr(0, i+1);
	prev.end_a = prev.start_a + prev.a.size() - q;
	prev.end_b = prev.start_b + prev.b.size() - r;
	prev.a = prev.a.substr(0, prev.a.size() - q);
	prev.b = prev.b.substr(0, prev.b.size() - r);
	// X="";for (auto c: prev.align_a) if (c!='-')X+=c; assert(X==prev.a);
	// assert(prev.a==qstr.substr(prev.start_a, prev.end_a-prev.start_a));
	// X="";for (auto c: prev.align_b) if (c!='-')X+=c; assert(X==prev.b);
	// assert(prev.b==rstr.substr(prev.start_b, prev.end_b-prev.start_b));



	q = 0, r = 0, i;
	for (i = 0; i < cur.alignment.size() && r < trim; i++) {
		if (cur.align_a[i] != '-') q++;
		if (cur.align_b[i] != '-') r++;
	}
	cur.align_a = cur.align_a.substr(i);
	cur.alignment = cur.alignment.substr(i);
	cur.align_b = cur.align_b.substr(i);
	cur.start_a += q;
	cur.start_b += r;
	cur.a = cur.a.substr(q);
	cur.b = cur.b.substr(r);
	// X="";for (auto c: cur.align_a) if (c!='-')X+=c; assert(X==cur.a);
		// assert(cur.a==qstr.substr(cur.start_a, cur.end_a-cur.start_a));
	// X="";for (auto c: cur.align_b) if (c!='-')X+=c; assert(X==cur.b);
		// assert(cur.b==rstr.substr(cur.start_b, cur.end_b-cur.start_b));
	
	prev.cigar_from_alignment();
	// eprn("")
	// eprn
	// eprn("{}\n{}", prev.align_a, prev.align_b);
		// eprn("{}", prev.print(70));

	cur.cigar_from_alignment();
		// eprn("{}", cur.print(70));


	// eprn("postfix... {}..{} to {}..{}", prev.start_a, prev.end_a, prev.start_b, prev.end_b);
	// eprn("         . {}..{} to {}..{}", cur.start_a, cur.end_a, cur.start_b, cur.end_b);

	assert(prev.start_a <= cur.start_a);
	assert(prev.start_b <= cur.start_b);
	assert(prev.end_a <= cur.start_a);
	assert(prev.end_b <= cur.start_b);
	int qgap = cur.start_a - prev.end_a;
	int rgap = cur.start_b - prev.end_b;
	// eprn("aligning {}..{} to {}..{}", prev.end_a, prev.end_a+qgap, prev.end_b, prev.end_b+rgap);
	// eprn("aligning {}..{} to {}..{}", prev.end_a, prev.end_a+qgap, prev.end_b, prev.end_b+rgap);
	if (qgap && rgap) {
		// eprn("aligning len {} -> {}", qgap, rgap);
		if (qgap <= 1000 && rgap <= 1000) { // "close" hits
			auto gap = align(qstr.substr(prev.end_a, qgap), rstr.substr(prev.end_b, rgap), 5, -4, 40, 1);
			prev.append_cigar(gap.cigar);
		} else { // assume only one part is the gap
			int ma = max(qgap, rgap); // gap 
			int mi = min(qgap, rgap); // mismatch
			// eprn("2xaligning len {}", mi);
			auto ma1 = align(qstr.substr(prev.end_a, mi), rstr.substr(prev.end_b, mi), 5, -4, 40, 1);
			ma1.cigar.push_back({qgap == mi ? 'I' : 'D', ma - mi});
			auto ma2 = align(qstr.substr(cur.start_a - mi, mi), rstr.substr(cur.start_b - mi, mi), 5, -4, 40, 1);
			ma2.cigar.push_front({qgap == mi ? 'I' : 'D', ma - mi});
			prev.append_cigar(ma2.error.error() < ma2.error.error() ? ma2.cigar : ma1.cigar);
		}
	} else if (qgap) {
		prev.append_cigar({{'D', qgap}});
	} else if (rgap) {
		prev.append_cigar({{'I', rgap}});	
	} 


	// eprn("done");
	prev.a += qstr.substr(prev.end_a, qgap) + cur.a; 
	prev.b += rstr.substr(prev.end_b, rgap) + cur.b; 
	assert(cur.end_a >= prev.end_a);
	assert(cur.end_b >= prev.end_b);
	prev.end_a = cur.end_a;
	prev.end_b = cur.end_b;
	prev.append_cigar(cur.cigar);
	prev.populate_nice_alignment();
	// eprn("{}", prev.print(70));
	// eprn("became... {}..{} to {}..{} --> {}", prev.start_a, prev.end_a, prev.start_b, prev.end_b, prev.cigar_string());

	// X="";for (auto c: prev.align_a) if (c!='-')X+=c; assert(X==prev.a); 
	// assert(prev.a==qstr.substr(prev.start_a, prev.end_a-prev.start_a));
	// X="";for (auto c: prev.align_b) if (c!='-')X+=c; assert(X==prev.b);
	// assert(prev.b==rstr.substr(prev.start_b, prev.end_b-prev.start_b));
}


void refine_chains(vector<Hit> &anchors, const string &qseq, const string &rseq)
{
	const double MATCH = 10;
	const double MISMATCH = 1;
	const double GAP = 0.5;
	const double GAPOPEN = 100;
	const double SIZE = 0;

	/// OVERLAPS!
	dprn(":: taking {} anchors for refinement", anchors.size());

	sort(anchors.begin(), anchors.end());

	vector<int> score;	
	for (auto &a: anchors) {
		score.push_back(
			+ MATCH * a.aln.error.matches 
			+ SIZE * a.aln.alignment.size()
			- MISMATCH * a.aln.error.mismatches 
			// - GAPOPEN * a.aln.error.gaps 
			- GAP * a.aln.error.gap_bases
		);
		dprn("-- init {}: (len:{}) {}..{} --> {}..{} ; score = {} cigar = {}", 
			&a-&anchors[0],
			abs(a.query_start- a.query_end), 
			a.query_start, a.query_end, a.ref_start, a.ref_end, score.back(),
			a.aln.cigar_string());
	}

	const int max_space = 10000;
	vector<int> dp(anchors.size(), 0);
	// vector<int> gaps_so_far(anchors.size(), 0);
	vector<int> prev(anchors.size(), -1);
	set<pair<int, int>, greater<pair<int, int>>> maxes;
	for (int ai = 0; ai < anchors.size(); ai++) {
		dp[ai] = score[ai];
		for (int aj = ai - 1; aj >= 0; aj--) {
			auto &c = anchors[ai];
			auto &p = anchors[aj];
			
			int cqs = c.query_start;
			if (cqs < p.query_end)
				cqs = p.query_end;
			int crs = c.ref_start;
			if (crs < p.ref_end)
				crs = p.ref_end;	

			if (p.query_end >= c.query_end || p.ref_end >= c.ref_end)
				continue;
			if (p.ref_start >= c.ref_start)
				continue;

			int ma = max(cqs - p.query_end, crs - p.ref_end);
			int mi = min(cqs - p.query_end, crs - p.ref_end);

			if (ma >= max_space)
				continue;

			int mis = MISMATCH * mi, //, GAPOPEN + GAP * 2 * mi),
				gap = GAPOPEN + GAP * (ma - mi);
			int sco = dp[aj] + score[ai] - mis - gap + SIZE * ma;
			// eprn("{} {} -> {} {} via {}", ai, aj, dp[ai], dp[aj], sco);
			if (sco >= dp[ai]) {
				dp[ai] = sco;
				prev[ai] = aj;
			}
		}
		maxes.insert({dp[ai], ai});
	}

	vector<bool> used(anchors.size(), 0);
	vector<deque<int>> paths;
	vector<Hit> hits;
	for (auto &m: maxes) {
		if (m.first == 0) 
			break;
		int maxi = m.second;
		if (used[maxi])
			continue;
		paths.push_back(deque<int>());
		while (maxi != -1 && (!used[maxi])) {
			paths.back().push_front(maxi);
			used[maxi] = true;
			maxi = prev[maxi];
		}

		int qlo = anchors[paths.back().front()].query_start, 
			qhi = anchors[paths.back().back()].query_end;
		int rlo = anchors[paths.back().front()].ref_start,   
			rhi = anchors[paths.back().back()].ref_end;

		// if (min(qhi - qlo, rhi - rlo) < MIN_READ_SIZE)
			// continue;

		int est_size = anchors[paths.back()[0]].aln.alignment.size();
		for (int i = 1; i < paths.back().size(); i++) {
			est_size += anchors[paths.back()[i]].aln.alignment.size();
			est_size += 
			max(anchors[paths.back()[i]].query_start - anchors[paths.back()[i-1]].query_end,
				anchors[paths.back()[i]].ref_start - anchors[paths.back()[i-1]].ref_end);
		}
		dprn("-- chain: (est:{}, size:{}) {}..{} --> {}..{} dp={}", est_size, 
			paths.back().size(), 
			qlo, qhi, rlo, rhi, dp[m.second]);
		for (auto p: paths.back()) {
			auto &y = anchors[p];
			dprn("    {}..{}->{}..{}", y.query_start, y.query_end, y.ref_start, y.ref_end);
		}
		auto hit = Hit {
			anchors.front().query, qlo, qhi,
			anchors.front().ref, rlo, rhi
		};


		if (est_size < 900)
			continue;
		vector<Hit> guide;
		Hit *prev = &anchors[paths.back()[0]];
		for (int pi = 1; pi < paths.back().size(); pi++) {
			auto &cur = anchors[paths.back()[pi]];
			if (cur.query_start < prev->query_end || cur.ref_start < prev->ref_end) {
				merge_alignments(prev->aln, cur.aln, qseq, rseq);

				// eprn("merging... {} to {}", &cur-&anchors[0], prev-&anchors[0]);
				// assert(prev->query_end <= cur.query_end);
				// assert(prev->ref_end <= cur.ref_end);
				
				// prev->query_end = cur.query_end;
				// prev->ref_end = cur.ref_end;
				// prev->aln = align(
				// 	hit.query->seq.substr(prev->query_start, 
				// 		prev->query_end - prev->query_start), 
				// 	hit.ref->seq.substr(prev->ref_start,
				// 		prev->ref_end - prev->ref_start),
				// 	5, -4, 40, 1);
				prev->query_start = prev->aln.start_a;
				prev->query_end = prev->aln.end_a;
				prev->ref_start = prev->aln.start_b;
				prev->ref_end = prev->aln.end_b;
				// eprn("{}", prev->aln.cigar_string());
			} else {
				guide.push_back(*prev);
				prev = &cur;
			}
		}
		guide.push_back(*prev);
		// eprn("done");

		hit.aln = Alignment::from_anchors(
			hit.query->seq, hit.ref->seq, guide, 250);
		// eprn("aln: {} {}", hit.aln.cigar_string(), hit.aln.alignment.size());
		hit.query_start = hit.aln.start_a; hit.ref_start = hit.aln.start_b;
		hit.query_end = hit.aln.end_a; hit.ref_end = hit.aln.end_b;
		// exit(0);
		if (hit.aln.alignment.size() >= 900)
			hits.push_back(hit);
	}

	anchors = hits;
}

// 310M1D393M1D27M2D170M4I1653M1D84M19D18M18D131M1D2542M1I33M43D4M43I364M111D8M111I39M5D2540M1D599M1D8204M153I121M78I47M232I28M78I230M155I307M1I538M1I153M1I262M1D416M5D104M4I70M2I741M1D510M1D357M1D63M4I483M1I945M1D415M2I2165M2I1835M1D843M1D411M1D310M1D191M1D3162M1D1485M4I1717M1I837M1D53M2D288M4I1427M2I1201M1D331M1D987M1I148M1I1514M190D190I1072M11I1199M1D748M3I632M5D53M6I2475M2D913M673D673I501M1D398M1D422M2I1345M133D1M133I13M2I1492M1D161M149D1M149I859M1I235M1D1466M53D53I1025M76D76I1012M44D44I97M1I975M3I463M23I74M3I200M2D1939M1I155M3I42M25D70M4I1337M2I808M1D141M5I156M29D28I482M1D144M4D25M1I350M4I268M1D1064M5D263M5D695M1I390M2I607M7D4587M1I38M4I388M2185I317M2694I337M2I210M2D282M1I27M


