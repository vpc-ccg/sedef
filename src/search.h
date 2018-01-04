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

typedef boost::icl::discrete_interval<int> Interval;
typedef boost::icl::interval_map<int, set<pair<Interval, Interval>>> Subtree;
typedef boost::icl::interval_map<int, Subtree> Tree;

/******************************************************************************/

struct Hit {
	int query_start, query_end; // query range
	int ref_start, ref_end; // reference range
	int jaccard; // coordinates of seed matches
	string reason;
};

/******************************************************************************/

vector<Hit> search (int query_winnow_start, 
	const Index &ref_hash, 
	const Index &query_hash, 
	Tree &tree,
	bool same_genome = true,
	int init_len = MIN_READ_SIZE,
	bool allow_extend = true,
	bool report_fails = false);

