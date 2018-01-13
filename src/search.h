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
#include "align.h"

/******************************************************************************/

typedef boost::icl::discrete_interval<int> Interval;
typedef boost::icl::interval_map<int, set<pair<Interval, Interval>>> Subtree;
typedef boost::icl::interval_map<int, Subtree> Tree;

/******************************************************************************/

vector<Hit> search (int query_winnow_start, 
	const Index &ref_hash, 
	const Index &query_hash, 
	Tree &tree,
	bool same_genome = true,
	int init_len = MIN_READ_SIZE,
	bool allow_extend = true,
	bool report_fails = false);

