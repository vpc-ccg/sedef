/// 786

/******************************************************************************/

#include <map>

#include "common.h"
#include "sliding.h"
#include "fasta.h"

using namespace std;

/******************************************************************************/

SlidingMap::SlidingMap(): 
	query_size(0), intersection(0), limit(0)
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

bool SlidingMap::add(const Hash &h, int BIT, int FULL) // bit: var with only one bit set
{
	auto it = this->lower_bound(h);
	
	bool inserted = false;
	if (it != this->end() && it->first == h) { 
		if (it->second & BIT) return false; 
		it->second |= BIT;
	} else {
		it = this->insert({h, BIT}).first;
		inserted = true;
	}

	assert(it != this->end());
	assert(boundary != this->end() || !query_size);
	if (query_size && it->first < boundary->first) {
		intersection += (it->second == FULL);
		if (inserted) {
			intersection -= (boundary->second == FULL);
			assert(boundary != this->begin()); // |S| >= 1!
			boundary--;
		}
	}
	return true;
}

bool SlidingMap::remove(const Hash &h, int BIT, int FULL) // bit: var with only one bit set
{	
	auto it = this->lower_bound(h); 
	if (it->first != h || !(it->second & BIT)) 
		return false;

	assert(boundary != this->end() || !query_size);
	if (query_size && it->first <= boundary->first) {
		intersection -= (it->second == FULL);
		if (it->second == BIT) {
			boundary++;
			assert(boundary != this->end());
			intersection += (boundary->second == FULL);
		}
	}

	assert(it != this->end());
	assert(it != boundary);
	if (it->second == BIT) {
		this->erase(it);
	} else {
		it->second &= ~BIT;
	}
	return true;
}

void SlidingMap::add_to_query(const Hash &h) 
{	
	if (!add(h, 1)) 
		return;

	limit = relaxed_jaccard_estimate(++query_size);
	assert(boundary != this->end() || (query_size == 1));
	if (boundary == this->end()) {
		boundary = this->begin();
	} else {
		boundary++;
	}
	assert(boundary != this->end());
	intersection += (boundary->second == 3);
}

void SlidingMap::remove_from_query(const Hash &h) 
{
	if (!remove(h, 1)) 
		return;
	
	limit = relaxed_jaccard_estimate(--query_size);
	assert(boundary != this->begin() || !query_size);
	intersection -= (boundary->second == 3);
	if (boundary == this->begin()) {
		boundary = this->end();
	} else { 
		boundary--;
	}
}

void SlidingMap::add_to_reference(const Hash &h) 
{	
	if (h.status != Hash::Status::HAS_N) {
		add(h, 2);
	}
}

void SlidingMap::remove_from_reference(const Hash &h) 
{
	if (h.status != Hash::Status::HAS_N) {
		remove(h, 2);
	}
}
