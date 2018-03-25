/// 786

/******************************************************************************/

#include <map>

#include "common.h"
#include "sliding.h"
#include "fasta.h"

using namespace std;

/******************************************************************************/

SlidingMap::SlidingMap(int kmer_size): 
	query_size(0), intersection(0), limit(0), kmer_size(kmer_size)
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
		// eprn("1");
		// eprn(">1> sz: {}..{}", this->size(), distance(this->begin(),boundary));
			boundary--;
		// eprn("<1< sz: {}..{}", this->size(), distance(this->begin(),boundary));

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
		// eprn("2");
		// eprn(">2> sz: {}..{}", this->size(), distance(this->begin(),boundary));
			boundary++;
		// eprn("<2< sz: {}..{}", this->size(), distance(this->begin(),boundary));
			assert(boundary != this->end());
			intersection += (boundary->second == FULL);
		}
	}

	assert(it != this->end());
	if (it->second == BIT) {
		assert(it != boundary);
		this->erase(it);
	} else {
		it->second &= ~BIT;
	}
	return true;
}

//// AFTER COPY
SlidingMap SlidingMap::fromMap(const SlidingMap &m)
{
	int dist = distance<SlidingMap::const_iterator>(m.begin(), m.boundary);
	SlidingMap s(m.kmer_size);
	s=m;
	s.boundary = s.begin();
	advance(s.boundary, dist);
	assert(m.boundary != m.end());
	assert(m.boundary->first == s.boundary->first);
	assert(m.boundary->second == s.boundary->second);
	return s;
}

void SlidingMap::add_to_query(const Hash &h) 
{	
	if (!add(h, 1)) 
		return;

	limit = relaxed_jaccard_estimate(++query_size, kmer_size) ;
	assert(boundary != this->end() || (query_size == 1));
	if (boundary == this->end()) {
		boundary = this->begin();
	} else {
		// eprn("3");
		// eprn(">3> sz: {}..{}", this->size(), distance(this->begin(),boundary));
		assert(boundary != this->end());
		boundary++;
		// eprn("<3< sz: {}..{}", this->size(), distance(this->begin(),boundary));

		// eprn("wooo");
	}
	assert(boundary != this->end());
	intersection += (boundary->second == 3);
}

void SlidingMap::remove_from_query(const Hash &h) 
{
	if (!remove(h, 1)) 
		return;
	
	limit = relaxed_jaccard_estimate(--query_size, kmer_size) ;
	assert(boundary != this->begin() || !query_size);
	intersection -= (boundary->second == 3);
	if (boundary == this->begin()) {
		boundary = this->end();
	} else { 
		// eprn("4");
		// eprn(">4> sz: {}..{}", this->size(), distance(this->begin(),boundary));
		boundary--;
		// eprn("<4< sz: {}..{}", this->size(), distance(this->begin(),boundary));
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
