/// 786

/******************************************************************************/

#include <map>

#include "common.h"
#include "sliding.h"

using namespace std;

/******************************************************************************/

SlidingMap::SlidingMap(): 
	query_size(0), time(0), intersection(0), limit(0) 
{ 
	boundary = this->end();
} 

int SlidingMap::jaccard() 
{
	if (intersection >= limit) {
		return intersection;
	} else {
		return intersection - limit;
	}
}

void SlidingMap::rewind()
{
	boundary = this->begin();
	std::advance(boundary, query_size - 1);
	assert(boundary != this->end());
	intersection = (boundary->second == 3);
	for (auto it = this->begin(); it != boundary; it++)
		intersection += (it->second == 3);
}

/// &1 = query; &2 = reference

void SlidingMap::add_to_query(const hash_t &h) 
{
	query_size++;
	limit = relaxed_jaccard_estimate(query_size);

	if (boundary != this->end() && next(boundary) != this->end()) {
		boundary++;
		intersection += (boundary->second == 3);
	}
	// if (h.first) return;

	auto it = this->lower_bound({h, 0}); 
	for (; it != this->end() && it->first.first == h && (it->second & 1); it++);

	bool inserted = false;
	if (it != this->end() && it->first.first == h && !(it->second & 1)) {
		it->second |= 1;
	} else {
		it = this->insert({{h, time++}, 1}).first;
		inserted = true;
	}

	if (boundary != this->end() && it->first <= boundary->first) {
		intersection += (it->second == 3);
		if (inserted) {
			intersection -= (boundary->second == 3);
			boundary--;
		}
	}
}

void SlidingMap::remove_from_query(const hash_t &h) // removing last added h from query!
{
	query_size--;
	limit = relaxed_jaccard_estimate(query_size);

	auto it = this->lower_bound({h, time + 1}); // first >=, and not = since...
	assert(it != this->begin());
	it--;
	while (it != this->begin() && it->first.first == h && !(it->second & 1)) it--; // find first recently added
	assert(it->first.first == h && (it->second & 1));
	if (it->second == 1) this->erase(it);
	else it->second &= 2;

	rewind();
}

void SlidingMap::remove_from_reference(const hash_t &h) // removing last added h from query!
{
	auto it = this->lower_bound({h, time + 1}); // first >=, and not = since...
	assert(it != this->begin());
	it--;
	while (it != this->begin() && it->first.first == h && !(it->second & 2)) it--; // find first recently added
	assert(it->first.first == h && (it->second & 2));
	if (it->second == 2) this->erase(it);
	else it->second &= 1;

	rewind();
}



void SlidingMap::add_to_reference(const hash_t &h)
{        
	if (h.first) return;
	
	// Returns an iterator pointing to the first element that is not less than key
	auto it = this->lower_bound({h, 0}); 
	for (; it != this->end() && it->first.first == h && (it->second & 2); it++);
	//     eprnn("// {}:{}:{} ", it->first.first.first, it->first.first.second, it->first.second);
	// }eprn("");

	bool inserted = false;
	if (it != this->end() && it->first.first == h && !(it->second & 2)) {
		it->second |= 2;
	} else {
		it = this->insert({make_pair(h, time++), 2}).first;
		inserted = true;
	}

	if (boundary != this->end() && it->first <= boundary->first) { // error
		intersection += (it->second == 3);
		if (inserted) {
			intersection -= (boundary->second == 3);
			boundary--;
		}
	}
}

void SlidingMap::remove_oldest_from_reference(const hash_t &h)
{
	if (h.first) return;
	
	auto it = this->lower_bound({h, 0}); 
	for (; it != this->end() && it->first.first == h && !(it->second & 2); it++);
	if (it == this->end() || it->first.first != h || !(it->second & 2)) return;

	it->second &= 1;
	if (boundary != this->end() && it->first <= boundary->first) {
		intersection -= (it->second == 1);
		if (it->second == 0) { // remove
			assert(next(boundary) != this->end());
			if (it != boundary) { /* make sure it is not boundary!! */
				this->erase(it);
				boundary++;
			} else {
				boundary++;
				this->erase(it);
			}
			intersection += (boundary->second == 3);
		}
	} else if (it->second == 0) {
		this->erase(it);
	}
}

