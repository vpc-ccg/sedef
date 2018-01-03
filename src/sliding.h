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

// struct _X {
// 	int
// };

struct SlidingMap: public std::map<hash_t, char> {
	typename SlidingMap::iterator boundary; // this is inclusive!
	int query_size;
	int intersection;
	double limit;
	int64_t time;

	SlidingMap();

	int jaccard();

	// when adding, if same hash is found: try to match earliest one if added by another set (i.e. try increase jaccard)
	// when removing, if same hash is found: try to remove the latest one (try to preserve jaccard)

	bool add(const hash_t &h, int BIT, int FULL = 3);
	bool remove(const hash_t &h, int BIT, int FULL = 3);

	void add_to_query(const hash_t &h);
	void remove_from_query(const hash_t &h);
	void add_to_reference(const hash_t &h);
	void remove_from_reference(const hash_t &h);
};
