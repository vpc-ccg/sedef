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

struct SlidingMap: public std::map<Hash, char> {
	typename SlidingMap::iterator boundary; // this is inclusive!
	int query_size;
	int intersection;
	double limit;
	int kmer_size;

	SlidingMap(int kmer_size);

	int jaccard();
	
	static SlidingMap fromMap(const SlidingMap &m);

	// when adding, if same hash is found: try to match earliest one if added by another set (i.e. try increase jaccard)
	// when removing, if same hash is found: try to remove the latest one (try to preserve jaccard)

	bool add(const Hash &h, int BIT, int FULL = 3);
	bool remove(const Hash &h, int BIT, int FULL = 3);

	void add_to_query(const Hash &h);
	void remove_from_query(const Hash &h);
	void add_to_reference(const Hash &h);
	void remove_from_reference(const Hash &h);
};
