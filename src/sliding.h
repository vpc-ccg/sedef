/// 786

/******************************************************************************/

#pragma once

/******************************************************************************/

#include <map>
#include <set>
#include <vector>
#include <string>
#include <queue>

#include "hash.h"

/******************************************************************************/

struct SlidingMap: public std::map<std::pair<hash_t, int64_t>, char> {
	typename std::map<std::pair<hash_t, int64_t>, char>::iterator boundary;
	int query_size;
	int intersection;
	double limit;
	int64_t time;

	SlidingMap();

	int jaccard();
	void rewind();

	// when adding, if same hash is found: try to match earliest one if added by another set (i.e. try increase jaccard)
	// when removing, if same hash is found: try to remove the latest one (try to preserve jaccard)

	void add_to_query(const hash_t &h);
	void remove_from_query(const hash_t &h);
	void add_to_reference(const hash_t &h);
	void remove_from_reference(const hash_t &h);
};

