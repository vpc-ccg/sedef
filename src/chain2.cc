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

void refine_chains(vector<Hit> &anchors)
{
	const double MATCH = 10;
	const double MISMATCH = 1;
	const double GAP = 1;
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
			- GAPOPEN * a.aln.error.gaps 
			- GAP * a.aln.error.gap_bases
		);
		// dprn("-- init {}: (len:{}) {}..{} --> {}..{} ; score = {} cigar = {}", 
		// 	&a-&anchors[0],
		// 	abs(a.query_start- a.query_end), 
		// 	a.query_start, a.query_end, a.ref_start, a.ref_end, score.back(),
		// 	a.aln.cigar_string());
	}

	const int max_space = 10000;
	vector<int> dp(anchors.size(), 0);
	vector<int> gaps_so_far(anchors.size(), 0);
	vector<int> prev(anchors.size(), -1);
	set<pair<int, int>, greater<pair<int, int>>> maxes;
	for (int ai = 0; ai < anchors.size(); ai++) {
		dp[ai] = score[ai];
		for (int aj = ai - 1; aj >= 0; aj--) {
			auto &c = anchors[ai];
			auto &p = anchors[aj];
			
			int cqs = c.query_start;
			// if (cqs < p.query_end)
			// 	cqs = p.query_end;
			int crs = c.ref_start;
			// if (crs < p.ref_end)
			// 	crs = p.ref_end;	

			if (cqs < p.query_end || crs < p.ref_end)
				continue;

			int ma = max(cqs - p.query_end, crs - p.ref_end);
			int mi = min(cqs - p.query_end, crs - p.ref_end);

			if (ma >= max_space)
				continue;

			int mis = min(MISMATCH * mi, GAPOPEN + GAP * 2 * mi),
				gap = GAPOPEN + GAP * (ma - mi);
			int sco = dp[aj] + score[ai] - mis - gap + SIZE * ma;
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
		// dprn("-- chain: (len:{}, size:{}, gaps={}) {}..{} --> {}..{} dp={}", abs(qhi-qlo), 
		// 	paths.back().size(), gaps,
		// 	qlo, qhi, rlo, rhi, dp[m.second]);
		auto hit = Hit {
			anchors.front().query, qlo, qhi,
			anchors.front().ref, rlo, rhi,
			0, "", "", {}, {}
		};

		if (est_size < 900)
			continue;

		for (auto p: paths.back()) {
			hit.guides.first.insert(hit.guides.first.end(), 
				anchors[p].guides.first.begin(), anchors[p].guides.first.end());
			hit.guides.second.insert(hit.guides.second.end(), 
				anchors[p].guides.second.begin(), anchors[p].guides.second.end());
		}

		hit.aln = Alignment::from_anchors(
			hit.query->seq, hit.ref->seq, hit.guides.first, hit.guides.second, 0);
		hit.query_start = hit.aln.start_a;
		hit.ref_start = hit.aln.start_b;
		hit.query_end = hit.aln.end_a;
		hit.ref_end = hit.aln.end_b;

		if (hit.aln.alignment.size() >= 900)
			hits.push_back(hit);
	}

	anchors = hits;
}

// Looking for 18,024..19,126 -> 17,524..19,268
//   Q -> GCCACCACCATGCCCAAGAG...CTGGATACAGTTTGGCTTGG
//   R -> GCCACCACCATGCCCAAGGG...ctggatacagttaggcttgg
//  18024.. 18394~>34=1   18394.. 18422~> 3=0   18422.. 18799~>34=1   18799.. 18808~> 1=0   18808.. 19118~>28=1   19118.. 19126~> 1=0  ; === 97
//  17524.. 17892~>21=1   17892.. 18232~>19=0   18232.. 18612~>22=1   18612.. 18950~>19=0   18950.. 19260~>18=1   19260.. 19268~> 0=0  ; === 97
