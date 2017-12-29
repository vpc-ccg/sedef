/// 786

/******************************************************************************/

#pragma once

/******************************************************************************/

#include <map>
#include <set>
#include <vector>
#include <string>
#include <queue>

#include <boost/icl/interval_map.hpp>
#include <boost/icl/interval_set.hpp>

#include "hash.h"

/******************************************************************************/

typedef boost::icl::discrete_interval<int> INTERVAL;
typedef boost::icl::interval_map<int, set<pair<INTERVAL, INTERVAL>>> SUBTREE_t;
typedef boost::icl::interval_map<int, SUBTREE_t> TREE_t;

/******************************************************************************/

struct Hit {
	int p, q; // query range
	int i, j; // reference range
	int jaccard; // coordinates of seed matches
	string reason;
};

/******************************************************************************/

vector<Hit> search (int query_start, 
	const Hash &ref_hash, 
	const Hash &query_hash, 
	TREE_t &tree,
	bool allow_overlaps = false,
	int init_len = MIN_READ_SIZE,
	bool allow_extend = true,
	bool report_fails = false);

