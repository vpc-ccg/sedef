/// 786
/// Adapted from SCALCE source code 
/// https://github.com/sfu-compbio/scalce

/******************************************************************************/

#pragma once

/******************************************************************************/

#include <list>
#include <vector>
#include <string>
#include <cmath>

#include "common.h"
#include "hit.h"

using namespace std;

/******************************************************************************/

struct Anchor {
	int query_start, query_end;
	int ref_start, ref_end;
	list<pair<int, int>> query_kmers, ref_kmers;

	bool operator< (const Anchor &a) const {
		return tie(query_start, ref_start) < tie(a.query_start, a.ref_start);
	}
};

std::vector<Hit> fast_align(const std::string &query, const std::string &ref);
