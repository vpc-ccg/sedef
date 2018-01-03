/// 786
/// Based on http://www.biorxiv.org/content/biorxiv/early/2017/03/24/103812.full.pdf

/******************************************************************************/

#include <list>
#include <queue>
#include <vector>
#include <string>
#include <cmath>
#include <queue>
#include <chrono>

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

const int MAX_MATCH = 1000 * 1000;   /// 1MB at most

/******************************************************************************/

#define eprn(f, ...)   // fmt::print(stderr, f "\n",  ##__VA_ARGS__)
#define eprnn(...)     // fmt::print(stderr, __VA_ARGS__)

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

bool is_in_tree(TREE_t &tree, TREE_t::const_iterator &pf, int pf_pos, int pfp_pos) 
{
	SUBTREE_t::const_iterator pfp;
	return (pf != tree.end() && (pfp = pf->second.find(pfp_pos)) != pf->second.end());
}

auto parse_hits(vector<Hit> &hits)
{
	vector<Hit> hits_real;
	for (auto &h: hits) { // TODO fix brute force
		bool add = true;
		for (auto &ph: hits) if ((&h - &hits[0]) != (&ph - &hits[0])) {
			if (h.i >= ph.i && h.j <= ph.j && h.p >= ph.p && h.q <= ph.q) { // if full match
				add = false;
				break;
			}
		}
		if (add) hits_real.push_back(h);
	}
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
	TREE_t &tree,
	bool report_fails) 
{
	eprn(">> extend query:{}..{} vs ref:{}..{}", query_start, query_end, ref_start, ref_end);
	
	assert(query_start < query_hash.seq.size());
	assert(ref_start < ref_hash.seq.size());
	assert(query_end <= query_hash.seq.size());
	assert(ref_end <= ref_hash.seq.size());

	auto f = filter(query_hash.seq, query_start, query_end - query_start, ref_hash.seq, ref_start, ref_end - ref_start);
	if (!f.first) {
		if (report_fails) hits.push_back({
			query_start, query_end, 
			ref_start, ref_end, 
			0,
			f.second
		});
		eprn(">> failed first filter {}", f.second);
		return;
	}

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

	// TODO: speed up by moving directly to next winnow

	for (int i = 0, j = winnow.jaccard(); ;) {
		int max_match = MAX_MATCH;
		if (!allow_overlaps) {
			max_match = min(max_match, int((1.0 / MAX_GAP_ERROR + .5) * abs(query_start - ref_start)));
		}
		if (max(query_end - query_start, ref_end - ref_start) > max_match) {
			break;
		}
		if (min(query_end - query_start, ref_end - ref_start) / (double)max(query_end - query_start, ref_end - ref_start) < (1 - 2 * MAX_GAP_ERROR)) {
			break;
		}

		bool succeeded = false;
		for (auto &fn: fns) {
			if ((i = fn.first()) == 0) 
				continue;
			if ((j = winnow.jaccard()) >= 0) {
				succeeded = true;
				break;
			} else {
				fn.second(i); // undo
			}
		}
		if (!succeeded) break;
	}
				
	eprn(">> extended to {}..{}; {}..{}", query_start, query_end, ref_start, ref_end);
	eprn(">>     took {}s", chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - time).count() / 1000.00);

	f = filter(query_hash.seq, query_start, query_end - query_start, ref_hash.seq, ref_start, ref_end - ref_start);
	if (!f.first) {
		if (report_fails) hits.push_back({
			query_start, query_end, 
			ref_start, ref_end, 
			winnow.jaccard(),
			f.second
		});
		eprn(">> failed second filter {}", f.second);
		return;
	}
	
	hits.push_back({
		query_start, query_end, 
		ref_start, ref_end, 
		winnow.jaccard(),
		"OK"
	});

	eprn(">> success!");
	
	auto a = INTERVAL(query_start, query_end);
	auto b = INTERVAL(ref_start, ref_end);
	tree += make_pair(a, SUBTREE_t({b, {make_pair(a, b)}}));
}

/******************************************************************************/

vector<Hit> search (int query_start, 
					const Hash &ref_hash, 
					const Hash &query_hash, 
					TREE_t &tree, 
					bool allow_overlaps,
					int init_len,
					bool allow_extend,
					bool report_fails)
{ 
	// TREE gives overlaps
	//  query vs ref

	if (query_start + init_len > query_hash.seq.size())
		return vector<Hit>();

	// if (memorized.size()) {
	// for (auto &me: memorized) {
	// 	prnn("{}-{} ", me.first, me.second);
	// } prn("");
	// }

	// CAN BE OPTIMIZED
	int st = query_hash.find_minimizers(query_start);
	if (st == query_hash.minimizers.size()) 
		return vector<Hit>();
	int mi = st;
	
	SlidingMap init_winnow;
	set<int> candidates_prel;
	auto pf = tree.find(query_start);
	for (; mi < query_hash.minimizers.size() && query_hash.minimizers[mi].second - query_start <= init_len; mi++) { 
		auto &h = query_hash.minimizers[mi].first;
		init_winnow.add_to_query(h);
		if (h.first != 0) // use only hashes with uppercase character!
			continue; 
		
		auto ptr = ref_hash.index.find(h);
		if (ptr == ref_hash.index.end() || ptr->second.size() >= ref_hash.threshold) {
			continue;
		} else for (auto &pos: ptr->second) {
			if (!allow_overlaps && pos < query_start + init_len) { // Make sure to have at least 1 kb spacing if reference = query
				continue;
			}
			if (!is_in_tree(tree, pf, query_hash.minimizers[mi].second, pos)) {
				candidates_prel.insert(pos);
			}
		}
	}
	if (!init_winnow.query_size)
		return vector<Hit>();
	int M = relaxed_jaccard_estimate(init_winnow.query_size);

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
		for (; ref_winnow_end < ref_hash.minimizers.size() && ref_hash.minimizers[ref_winnow_end].second < ref_end; ref_winnow_end++) 
			winnow.add_to_reference(ref_hash.minimizers[ref_winnow_end].first);
		// winnow.rewind();

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
			if (ref_winnow_start < ref_hash.minimizers.size() && ref_hash.minimizers[ref_winnow_start].second < ref_start + 1)
				winnow.remove_from_reference(ref_hash.minimizers[ref_winnow_start++].first);
			if (ref_winnow_end < ref_hash.minimizers.size() && ref_hash.minimizers[ref_winnow_end].second == ref_end)
				winnow.add_to_reference(ref_hash.minimizers[ref_winnow_end++].first);
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
			if (report_fails) hits.push_back({
				query_start, query_start + init_len, 
				best_ref_start, best_ref_end, 
				best_jaccard, 
				fmt::format("jaccard: {} < {}", best_winnow.limit + best_jaccard, best_winnow.limit)
			});
			JACCARD_FAILED++;
		} else if (allow_extend) {
			if (!check_overlap(tree, pf, query_start, query_start + init_len, best_ref_start, best_ref_end)) {
				INTERVAL_FAILED++;
			} else {
				// eprn("init jacc === {}", best_jaccard);
				extend(best_winnow,
					query_hash, query_start, query_start + init_len, st, mi,
					ref_hash, best_ref_start, best_ref_end, best_ref_winnow_start, best_ref_winnow_end, 
					hits, allow_overlaps, tree, report_fails
				);
				// pf = tree.find(query_start);
			}
		} else {
			auto f = filter(query_hash.seq, query_start, init_len, ref_hash.seq, best_ref_start, best_ref_end - best_ref_start);
			if (f.first || report_fails) hits.push_back({
				query_start, query_start + init_len, 
				best_ref_start, best_ref_end, 
				best_jaccard, 
				f.second == "" ? "OK_INIT" : f.second
			});
		}
	}
   
	tree -= INTERVAL(0, query_start - MIN_READ_SIZE);
	return parse_hits(hits);
}
