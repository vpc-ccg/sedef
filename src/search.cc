/// 786
/// Based on http://www.biorxiv.org/content/biorxiv/early/2017/03/24/103812.full.pdf

/******************************************************************************/

#include <list>
#include <queue>
#include <vector>
#include <string>
#include <cmath>
#include <queue>

#include "common.h"
#include "filter.h"
#include "search.h"
#include "sliding.h"

using namespace std;

/******************************************************************************/

/* extern */ int64_t TOTAL_ATTEMPTED = 0;
/* extern */ int64_t JACCARD_FAILED = 0;
/* extern */ int64_t INTERVAL_FAILED = 0;

/******************************************************************************/

const int MAX_MATCH = 500 * 1000;   /// 50kb at most
const int RECOVER_BP = 1000;        /// We allow 250bp extra extend just in case!

/******************************************************************************/

bool check_overlap(TREE_t &tree, TREE_t::const_iterator &pf, int pf_pos, int pf_end, int pfp_pos, int pfp_end) 
{
	assert(pf_pos <= pf_end);
	assert(pfp_pos <= pfp_end);

	const int RIGHT_ALLOWANCE = 750;    /// TODO more mathy formulation

	SUBTREE_t::const_iterator pfp;
	if (pf != tree.end() && (pfp = pf->second.find(pfp_pos)) != pf->second.end()) { // check overlap!
		for (auto &it: pfp->second) {
			int sA = it.first.lower(),  eA = it.first.upper();
			int sB = it.second.lower(), eB = it.second.upper();

			if (pf_pos >= sA && pf_end <= eA && pfp_pos >= sB && pfp_end <= eB)
				return false;

			if (min(eA - sA, eB - sB) < MIN_READ_SIZE * 1.5) 
				continue;

			/*
			-------+     
			  +---------
			  <----> must be at least RIGHT_ALLOWANCE 
			*/
			if (eA - pf_pos  >= RIGHT_ALLOWANCE && eB - pfp_pos >= RIGHT_ALLOWANCE) 
				return false;
		}
	}
	return true;
}

auto parse_hits(vector<Hit> &hits)
{
	// COUNT SQUEEZED!
	vector<Hit> hits_real;
	// sort(hits.begin(), hits.end(), [](const Hit &a, const Hit &b) {
	//     // if (a.q != b.q) return a.q > b.q; // .query_start are all equal
	//     return make_pair(a.i, -a.j) < make_pair(b.i, -b.j);   // larger intervals at the top! then by other ints
	// });
	for (auto &h: hits) { // brute force, but maybe easier than doing full interval thing
		// eprn("{:>9} {:>9} --- {:>9} {:>9}", h.p, h.q, h.i, h.j);
		bool add = true;
		for (auto &ph: hits) if ((&h - &hits[0]) != (&ph - &hits[0])) {
			if (h.i >= ph.i && h.j <= ph.j && h.p >= ph.p && h.q <= ph.q) { // if full match
				add = false;
				break;
			}
		}
		if (add) hits_real.push_back(h);
	}
	//reverse(hits_real.begin(), hits_real.end()); // smaller to larger
	return hits_real;
}

/******************************************************************************/

void extend(SlidingMap &winnow,
	const Hash &query_hash, int query_start, int query_end,
	int query_winnow_start, int query_winnow_end,
	const Hash &ref_hash, int ref_start, int ref_end,
	int ref_winnow_start, int ref_winnow_end,
	vector<Hit> &hits,
	bool allow_overlaps,
	TREE_t &tree) 
{
	static int CNT(0);
	++CNT;

	eprn(">> {}: extend query:{}..{} vs ref:{}..{}", CNT, query_start, query_end, ref_start, ref_end);
	
	assert(query_start < query_hash.seq.size());
	assert(ref_start < ref_hash.seq.size());
	assert(query_end <= query_hash.seq.size());
	assert(ref_end <= ref_hash.seq.size());
	auto f = filter(query_hash.seq, query_start, query_end - query_start, ref_hash.seq, ref_start, ref_end - ref_start);
	if (!f.first) {
		hits.push_back({
			query_start, query_end, 
			ref_start, ref_end, 
			0,
			f.second
		});
		// eprn(">> {}: failed first filter {}", CNT, f.second);
		return;
	}

	int j = winnow.jaccard();
	auto best_winnow = winnow;
	int best_query_end = query_end, // TODO is not inclusive!
	    best_query_winnow_end = query_winnow_end,
	    best_query_start = query_start,
	    best_query_winnow_start = query_winnow_start,
	    best_ref_end = ref_end,
	    best_ref_winnow_end = ref_winnow_end,
	    best_ref_start = ref_start,
	    best_ref_winnow_start = ref_winnow_start,
	    best_jaccard = j;

	auto fn_qe = [&]() {
		if (query_end >= query_hash.seq.size()) return 0;
		int i = 1;
		if (query_winnow_end < query_hash.minimizers.size() && query_hash.minimizers[query_winnow_end].second == query_end)
			winnow.add_to_query(query_hash.minimizers[query_winnow_end++].first), i++;
		query_end++;
		return i;
	};
	auto fn_u_qe = [&](int i) {
		if (!i) return;
		if (i == 2) winnow.remove_from_query(query_hash.minimizers[--query_winnow_end].first);
		query_end--;
	};

	auto fn_re = [&]() {
		if (ref_end >= ref_hash.seq.size()) return 0;
		int i = 1;
		if (ref_winnow_end < ref_hash.minimizers.size() && ref_hash.minimizers[ref_winnow_end].second == ref_end) 
			winnow.add_to_reference(ref_hash.minimizers[ref_winnow_end++].first), i++;
		ref_end++;
		return i;
	};
	auto fn_u_re = [&](int i) {
		if (!i) return;
		if (i == 2) winnow.remove_from_reference(ref_hash.minimizers[--ref_winnow_end].first);
		ref_end--;
	};

	auto fn_qre = [&]() {
		if (query_end >= query_hash.seq.size() || ref_end >= ref_hash.seq.size()) return 0;
		return fn_qe() * 10 + fn_re();
	};
	auto fn_u_qre = [&](int i) {
		if (!i) return;
		fn_u_qe(i / 10), fn_u_re(i % 10);
	};

	auto fn_qs = [&]() {
		if (!query_start) return 0;
		int i = 1;
		if (query_winnow_start && query_hash.minimizers[query_winnow_start - 1].second == query_start - 1)
			winnow.add_to_query(query_hash.minimizers[--query_winnow_start].first), i++;
		query_start--;
		return i;
	};
	auto fn_u_qs = [&](int i) {
		if (!i) return;
		if (i == 2) winnow.remove_from_query(query_hash.minimizers[query_winnow_start++].first);
		query_start++;
	};

	auto fn_rs = [&]() {
		if (!ref_start) return 0;
		int i = 1;
		if (ref_winnow_start && ref_hash.minimizers[ref_winnow_start - 1].second == ref_start - 1) 
			winnow.add_to_reference(ref_hash.minimizers[--ref_winnow_start].first), i++; 
		ref_start--;
		return i;
	};
	auto fn_u_rs = [&](int i) {
		if (!i) return;
		if (i == 2) winnow.remove_from_reference(ref_hash.minimizers[ref_winnow_start++].first);
		ref_start++;
	};

	auto fn_qrs = [&]() {
		if (!query_start || !ref_start) return 0;
		return fn_qs() * 10 + fn_rs();
	};
	auto fn_u_qrs = [&](int i) {
		if (!i) return;
		fn_u_qs(i / 10), fn_u_rs(i % 10);
	};
	
	auto fns = vector<pair<function<int(void)>, function<void(int)>>>{
		make_pair(fn_qre, fn_u_qre),
		make_pair(fn_qe, fn_u_qe),
		make_pair(fn_re, fn_u_re),
		make_pair(fn_qrs, fn_u_qrs),
		make_pair(fn_qs, fn_u_qs),
		make_pair(fn_rs, fn_u_rs)
	};

	const int STEP = 100;
	for ( ; ; ) {
		int max_match = MAX_MATCH;
		if (!allow_overlaps)
			max_match = min(max_match, int((1.0 / MAX_GAP_ERROR + .5) * abs(query_start - ref_start)));
		if (max(query_end - query_start, ref_end - ref_start) > max_match)
			break;
		if (min(query_end - query_start, ref_end - ref_start) / (double)max(query_end - query_start, ref_end - ref_start) < (1 - 2 * MAX_GAP_ERROR))
			break;

		bool succeeded = false;
		for (auto &fn: fns) {
			// eprn("{} , {}--{}/{}--{}", CNT, query_start, query_end, ref_start, ref_end);
			int i = fn.first();
			if (i == 0) continue;
			//for (i = 0; i < STEP; i++) {
			//}
			//if (i == 0) continue;
			// we have "working" jacard
			if ((j = winnow.jaccard()) >= 0) {
				// best_jaccard = j;
				// best_query_start = query_start;
				// best_query_winnow_start = query_winnow_start;
				// best_query_end = query_end;
				// best_query_winnow_end = query_winnow_end;
				// best_ref_start = ref_start;
				// best_ref_winnow_start = ref_winnow_start;
				// best_ref_end = ref_end;
				// best_ref_winnow_end = ref_winnow_end;
				// best_winnow = winnow;

				succeeded = true;
				break;
			} else {
				fn.second(i); // undo
			}
			// else {
			// 	query_start = best_query_start;
			// 	query_winnow_start = best_query_winnow_start;
			// 	query_end = best_query_end;
			// 	query_winnow_end = best_query_winnow_end;

			// 	ref_start = best_ref_start;
			// 	ref_winnow_start = best_ref_winnow_start;
			// 	ref_end = best_ref_end;
			// 	ref_winnow_end = best_ref_winnow_end;

			// 	winnow = best_winnow;
			// }
		}
		if (!succeeded) break;
	}

	best_query_start = query_start;
	best_query_end = query_end;
	best_ref_start = ref_start;
	best_ref_end = ref_end;
				
	f = filter(query_hash.seq, best_query_start, best_query_end - best_query_start, 
		ref_hash.seq, best_ref_start, best_ref_end - best_ref_start);
	if (!f.first) {
		hits.push_back({
			best_query_start, best_query_end, 
			best_ref_start, best_ref_end, 
			best_jaccard,
			f.second
		});
		// eprn(">> {}: failed second filter {}", CNT, f.second);
		return;
	}
	
	hits.push_back({
		best_query_start, best_query_end, 
		best_ref_start, best_ref_end, 
		best_jaccard,
		"OK"
	});

	// eprn(">> {}: extended to {}..{}; {}..{}", CNT, best_query_start, best_query_end, 
	// 	best_ref_start, best_ref_end);
	
	auto a = INTERVAL(query_start, best_query_end);
	auto b = INTERVAL(ref_start, best_ref_end);
	tree += make_pair(a, SUBTREE_t({b, {make_pair(a, b)}}));
}

/******************************************************************************/

vector<Hit> search (int query_start, 
					const Hash &ref_hash, 
					const Hash &query_hash, 
					TREE_t &tree, 
					bool allow_overlaps,
					int init_len)
{ 
	if (query_start + init_len > query_hash.seq.size())
		return vector<Hit>();

	int st = query_hash.find_minimizers(query_start);
	if (st == query_hash.minimizers.size()) 
		return vector<Hit>();
	int mi = st;
	
	auto pf = tree.find(query_start);

	// Iterate through all unique hashes in query[query_start: query_start + init_len]
	// `candidates` is a list of positions in the reference which match the query hashes
	SlidingMap init_winnow;
	set<int> candidates_prel;
	for (; mi < query_hash.minimizers.size() && query_hash.minimizers[mi].second - query_start <= init_len; mi++) { 
		auto &h = query_hash.minimizers[mi].first;
		init_winnow.add_to_query(h);

		if (h.first) // If it is N hash, ignore it
			continue;
		
		auto ptr = ref_hash.index.find(h);
		if (ptr == ref_hash.index.end() || ptr->second.size() >= ref_hash.threshold) 
			continue;
		for (auto &pos: ptr->second) {
			// Make sure to have at least 1 kb spacing if reference = query
			if (!allow_overlaps && pos < query_start + init_len)
				continue;
			// if (!check_overlap(pf, query_start, pos)) {
			//     continue;
			// }
			candidates_prel.insert(pos);
		}
	}
	if (!init_winnow.query_size) { // TODO min sketch size limit for Ns
		return vector<Hit>();
	}
	int MO = ceil(init_winnow.query_size * tau());
	int M  = relaxed_jaccard_estimate(init_winnow.query_size);

	// eprn("==MO {} M {}", MO, M);
	// int M = init_winnow.size() * tau();
	
	// for (auto e: candidates) eprn("{:n}", e);
	// eprn(">> search: {} candidates to evaluate (M={}, init_winnow_query={}, tau={})", candidates_prel.size(), M, init_winnow.query_size, tau());

	// Find all locations in the `candidates` so that
	// init_len read covers at least M (= s * tau) hashes
	//sort(candidates.begin(), candidates.end());
	vector<int> candidates(candidates_prel.begin(), candidates_prel.end());
	vector<pair<int, int>> T;
	for (int i = 0; i <= (int)candidates.size() - M; i++) {
		int j = i + (M - 1);
		if (candidates[j] - candidates[i] < init_len) {
			int x = max(0, candidates[j] - init_len + 1), 
				y = candidates[i] + 1;
			if (T.size() && x < T.back().second) {
				T.back().second = max(T.back().second, y);
			} else {
				T.push_back({x, y});
			}
		}
	}
	vector<Hit> hits;
	for (auto &t: T) {
		if (!allow_overlaps)
			t.first = max(t.first, query_start + init_len);
		if (t.first > t.second) 
			continue;

		// eprn(">> BEGIN extend query:{:n} vs ref:{:n}...{:n}", query_start, t.first, t.second);
	
		TOTAL_ATTEMPTED++; 
		auto winnow = init_winnow;
	 
		int ref_start = t.first, 
			ref_end   = min(t.first + init_len, (int)ref_hash.seq.size());
		assert(ref_start >= 0);
		int ref_winnow_start = ref_hash.find_minimizers(ref_start);
		assert(ref_winnow_start < ref_hash.minimizers.size());
		assert(winnow.query_size > 0);

		// winnow is W(query) ; extend it to W(query) | W(ref) and mark elements in W(query) & W(ref)
		int ref_winnow_end = ref_winnow_start;
		for (; ref_winnow_end < ref_hash.minimizers.size() && ref_hash.minimizers[ref_winnow_end].second < ref_end; ref_winnow_end++) {
			winnow.add_to_reference(ref_hash.minimizers[ref_winnow_end].first);
		}
		winnow.rewind();

		// eprn(">> extend needed_jaccard:{} jaccard:{} init_start:{} init_jaccard={}", 
		//     tau * winnow_query.size(),
		//     jaccard, 
		//     seed_ref_start, seed_jaccard);
		// Roll until we find best inital match
		auto best_winnow = winnow;
		int best_jaccard = winnow.jaccard();
		int best_ref_start = ref_start, best_ref_end = ref_end;
		int best_ref_winnow_start = ref_winnow_start, best_ref_winnow_end = ref_winnow_end;
		while (ref_start < t.second && ref_end < ref_hash.seq.size()) {
			if (ref_winnow_start < ref_hash.minimizers.size() && ref_hash.minimizers[ref_winnow_start].second < ref_start + 1) {
				winnow.remove_oldest_from_reference(ref_hash.minimizers[ref_winnow_start].first);
				ref_winnow_start++;    
			}
			if (ref_winnow_end < ref_hash.minimizers.size() && ref_hash.minimizers[ref_winnow_end].second < ref_end + 1) {
				winnow.add_to_reference(ref_hash.minimizers[ref_winnow_end].first);
				ref_winnow_end++;
			}
			int j;
			if ((j = winnow.jaccard()) > best_jaccard) {
				best_jaccard = j;
				best_ref_start = ref_start;
				best_ref_end = ref_end;
				best_ref_winnow_start = ref_winnow_start;
				best_ref_winnow_end = ref_winnow_end;
				best_winnow = winnow;
			}
			ref_start++;
			ref_end++;
			if (ref_end == ref_hash.seq.size()) 
				break;
		}

		if (best_jaccard < 0) {
			hits.push_back({
				query_start, query_start + init_len, 
				best_ref_start, best_ref_end, 
				best_jaccard, 
				fmt::format("jaccard: {} < {}", best_winnow.limit + best_jaccard, best_winnow.limit)
			});
			JACCARD_FAILED++;
		} else if (init_len == MIN_READ_SIZE) {
			auto pf = tree.find(query_start);
			if (!check_overlap(tree, pf, query_start, query_start + init_len, best_ref_start, best_ref_end)) {
				INTERVAL_FAILED++;
			} else {
				// eprn("init jacc === {}", best_jaccard);
				extend(best_winnow,
					query_hash, query_start, query_start + init_len, st, mi,
					ref_hash, best_ref_start, best_ref_end, best_ref_winnow_start, best_ref_winnow_end, 
					hits, allow_overlaps, tree
				); 
			}
		} else {
			auto f = filter(query_hash.seq, query_start, init_len, ref_hash.seq, best_ref_start, best_ref_end - best_ref_start);
			hits.push_back({
				query_start, query_start + init_len, 
				best_ref_start, best_ref_end, 
				best_jaccard, 
				f.second == "" ? "OK_INIT" : f.second
			});
		}
	}
   
	tree -= INTERVAL(0, query_start);
	return parse_hits(hits);
}

