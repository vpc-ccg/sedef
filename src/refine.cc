/// 786

#include "common.h"
#include "fasta.h"
#include "align.h"
#include "search.h"
#include "segment.h"
#include "chain.h"
#include "hit.h"
#include "refine.h"

using namespace std; 

/******************************************************************************/

const double REFINE_MATCH = 10;
const double REFINE_MISMATCH = 1;
const double REFINE_GAP = 0.5;
const double REFINE_GAPOPEN = 100;

const int REFINE_MIN_READ = 900;
const int REFINE_SIDE_ALIGN = 500;
const int REFINE_MAX_GAP = 10000;

/******************************************************************************/

void refine_chains(vector<Hit> &anchors, const string &qseq, const string &rseq,
	const Hit &orig)
{
	dprn(":: taking {} anchors for refinement", anchors.size());

	sort(anchors.begin(), anchors.end());

	bool same_chr = orig.query->name == orig.ref->name && orig.query->is_rc == orig.ref->is_rc;
	vector<int> score;	
	for (auto &a: anchors) {
		score.push_back(
			+ REFINE_MATCH * a.aln.matches()
			- REFINE_MISMATCH * a.aln.mismatches()
			- REFINE_GAP * a.aln.gap_bases()
		);

		// dprn("-- init {}: (len:{}) {}..{} --> {}..{} ; score = {} cigar = {}", 
		// 	&a-&anchors[0],
		// 	abs(a.query_start- a.query_end), 
		// 	a.query_start, a.query_end, a.ref_start, a.ref_end, score.back(),
		// 	a.aln.cigar_string());
	}


	vector<int> dp(anchors.size(), 0);
	vector<int> prev(anchors.size(), -1);
	set<pair<int, int>, greater<pair<int, int>>> maxes;
	for (int ai = 0; ai < anchors.size(); ai++) {
		if (same_chr) {
			auto &c = anchors[ai];
			int qlo = c.query_start, qhi = c.query_end;
			int rlo = c.ref_start, rhi = c.ref_end;
			int qo = max(0, 
					min(orig.query_start + qhi, orig.ref_start + rhi) - 
				 	max(orig.query_start + qlo, orig.ref_start + rlo));
			if ((rhi-rlo)-qo < REFINE_SIDE_ALIGN && (qhi-qlo)-qo<REFINE_SIDE_ALIGN) { // no gap between 
				// eprn("REJECTED {}..{} --> {}..{}",
					// c.query_start+orig.query_start, c.query_end+orig.query_start,
					// c.ref_start+orig.ref_start, c.ref_end+orig.ref_start);
				continue;
			}
		}

		dp[ai] = score[ai];
		for (int aj = ai - 1; aj >= 0; aj--) {
			auto &c = anchors[ai];
			auto &p = anchors[aj];
			
			int cqs = c.query_start;
			if (cqs < p.query_end) {
				cqs = p.query_end;
			}
			int crs = c.ref_start;
			if (crs < p.ref_end) {
				crs = p.ref_end;	
			}

			if (p.query_end >= c.query_end || p.ref_end >= c.ref_end)
				continue;
			if (p.ref_start >= c.ref_start)
				continue;

			int ma = max(cqs - p.query_end, crs - p.ref_end);
			int mi = min(cqs - p.query_end, crs - p.ref_end);

			if (ma >= REFINE_MAX_GAP)
				continue;

			if (same_chr) {
				int qlo = p.query_end, qhi = cqs;
				int rlo = p.ref_end, rhi = crs;
				int qo = max(0, 
						min(orig.query_start + qhi, orig.ref_start + rhi) - 
					 	max(orig.query_start + qlo, orig.ref_start + rlo));
				if (qo >= 1) { // no gap between 
					continue;
				}
			}

			int mis = REFINE_MISMATCH * mi, 
				gap = REFINE_GAPOPEN + REFINE_GAP * (ma - mi);
			int sco = dp[aj] + score[ai] - mis - gap;
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
		int hasu = 0;
		while (maxi != -1 && (!used[maxi])) {
			paths.back().push_front(maxi);
			hasu += anchors[maxi].jaccard;
			used[maxi] = true;
			maxi = prev[maxi];
		}

		int qlo = anchors[paths.back().front()].query_start, 
			qhi = anchors[paths.back().back()].query_end;
		int rlo = anchors[paths.back().front()].ref_start,   
			rhi = anchors[paths.back().back()].ref_end;

		// if (min(qhi - qlo, rhi - rlo) < MIN_READ_SIZE)
			// continue;

		int est_size = anchors[paths.back()[0]].aln.span();
		for (int i = 1; i < paths.back().size(); i++) {
			est_size += anchors[paths.back()[i]].aln.span();
			est_size += 
			max(anchors[paths.back()[i]].query_start - anchors[paths.back()[i-1]].query_end,
				anchors[paths.back()[i]].ref_start - anchors[paths.back()[i-1]].ref_end);
		}
		dprn("-- chain: (est:{}, size:{}) {}..{} --> {}..{}  ## {}..{} --> {}..{} ", est_size, 
			paths.back().size(), 
			qlo, qhi, rlo, rhi, //dp[m.second],
			qlo+orig.query_start, qhi+orig.query_start, rlo+orig.ref_start, rhi+orig.ref_start);
		for (auto p: paths.back()) {
			auto &y = anchors[p];
			dprn("    {}..{}->{}..{}", y.query_start, y.query_end, y.ref_start, y.ref_end);
		}

		if (est_size < REFINE_MIN_READ - REFINE_SIDE_ALIGN) {
			dprn("est size failed");
			continue;
		}


		// if (same_chr) {
		// 	// assert(orig.query->is_rc == 0);
		// 	int qo = max(0, 
		// 				min(orig.query_start + qhi, orig.ref_start + rhi) - 
		// 			 	max(orig.query_start + qlo, orig.ref_start + rlo));
		// 	if (pct(qo, rhi - rlo) || pct(qo, qhi - qlo)) {
		// 		dprn("overlap in genome failed");

		// 		// attempt to fix it
		// 		if 

		// 		continue;
		// 	}
		// }
		bool overlap = 0;
		for (auto &h: hits) {
			int qo = max(0, min(qhi, h.query_end) - max(qlo, h.query_start));
			int ro = max(0, min(rhi, h.ref_end) - max(rlo, h.ref_start));

			if (qhi-qlo-qo < REFINE_SIDE_ALIGN && rhi-rlo-ro < REFINE_SIDE_ALIGN) {
				dprn("between overlap failed");
				overlap = 1;
				break;
			}		
		}
		if (overlap)
			continue;

		auto hit = Hit {
			anchors.front().query, qlo, qhi,
			anchors.front().ref, rlo, rhi
		};

		vector<Hit> guide;
		Hit *prev = &anchors[paths.back()[0]];
		for (int pi = 1; pi < paths.back().size(); pi++) {
			auto &cur = anchors[paths.back()[pi]];
			if (cur.query_start < prev->query_end || cur.ref_start < prev->ref_end) {
				prev->aln.merge(cur.aln, qseq, rseq);
				update_from_alignment(*prev);
			} else {
				guide.push_back(*prev);
				prev = &cur;
			}
		}
		guide.push_back(*prev);

		hit.aln = Alignment(hit.query->seq, hit.ref->seq, guide, REFINE_SIDE_ALIGN);
		update_from_alignment(hit);
		// eprn("aln: {} {}", hit.aln.cigar_string(), hit.aln.span());		
		if (hit.aln.span() >= REFINE_MIN_READ) {
			dprn("IN!");
			hits.push_back(hit);
		} else {
			dprn("failed final size ");
		}
	}

	anchors = hits;
}

// 310M1D393M1D27M2D170M4I1653M1D84M19D18M18D131M1D2542M1I33M43D4M43I364M111D8M111I39M5D2540M1D599M1D8204M153I121M78I47M232I28M78I230M155I307M1I538M1I153M1I262M1D416M5D104M4I70M2I741M1D510M1D357M1D63M4I483M1I945M1D415M2I2165M2I1835M1D843M1D411M1D310M1D191M1D3162M1D1485M4I1717M1I837M1D53M2D288M4I1427M2I1201M1D331M1D987M1I148M1I1514M190D190I1072M11I1199M1D748M3I632M5D53M6I2475M2D913M673D673I501M1D398M1D422M2I1345M133D1M133I13M2I1492M1D161M149D1M149I859M1I235M1D1466M53D53I1025M76D76I1012M44D44I97M1I975M3I463M23I74M3I200M2D1939M1I155M3I42M25D70M4I1337M2I808M1D141M5I156M29D28I482M1D144M4D25M1I350M4I268M1D1064M5D263M5D695M1I390M2I607M7D4587M1I38M4I388M2185I317M2694I337M2I210M2D282M1I27M

