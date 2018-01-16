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

#include "hit.h"
#include "hash.h"
#include "align.h"

/******************************************************************************/

typedef boost::icl::discrete_interval<int> Interval;
typedef boost::icl::interval_map<int, set<pair<Interval, Interval>>> Subtree;
typedef boost::icl::interval_map<int, Subtree> Tree;

/******************************************************************************/

vector<Hit> search (int query_winnow_start, 
	shared_ptr<Index> query_hash, 
	shared_ptr<Index> ref_hash, 
	Tree &tree,
	const bool same_genome = true,
	const int init_len = MIN_READ_SIZE,
	const bool allow_extend = true,
	const bool report_fails = false);

