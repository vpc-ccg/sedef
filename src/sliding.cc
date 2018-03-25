/// 786

/******************************************************************************/

#include <map>

#include "common.h"
#include "sliding.h"
#include "fasta.h"

using namespace std;

/******************************************************************************/

SlidingMap::SlidingMap(int kmer_size): 
	query_size(0), intersection(0), limit(0), kmer_size(kmer_size), storage()
{ 
	boundary = storage.end();
	// dprn("> Constructor:{:5} > boundary is at the end", storage.size());
}

SlidingMap::SlidingMap(const SlidingMap &other):
	query_size(other.query_size), intersection(other.intersection), 
	limit(other.limit), kmer_size(other.kmer_size),
	storage(other.storage)
{
	if (other.boundary == other.storage.end()) {
		boundary = storage.end();
	} else {
		// dprn("> Copy:{:5} > boundary {} in orig is at position {}", 
		// 	other.storage.size(), other.print_it(other.boundary), 
		// 	distance<map<Hash, char>::const_iterator>(other.storage.begin(), other.boundary));

		boundary = storage.begin();
		int dist = distance<map<Hash, char>::const_iterator>(other.storage.begin(), other.boundary);
		advance(boundary, dist);
		assert(boundary != storage.end());
		assert(boundary->first == other.boundary->first);
		assert(boundary->second == other.boundary->second);

		// dprn("> Copy:{:5} > boundary {} in new is at position {}", 
		// 	storage.size(), print_it(boundary), 
		// 	distance(storage.begin(), boundary));
	}
}

SlidingMap::SlidingMap(SlidingMap&& other): SlidingMap()
{
	swap(*this, other);
}

SlidingMap& SlidingMap::operator=(SlidingMap other)
{
	swap(*this, other);
	return *this;
}

/******************************************************************************/

int SlidingMap::jaccard() 
{
	if (intersection >= limit) {
		return intersection;
	} else {
		return intersection - limit;
	}
}

string SlidingMap::print_it(const map<Hash,char>::iterator &b) const
{
	if (b == storage.end()) {
		return "END";
	} else {
		return fmt::format("{}:{}", b->first, b->second);
	}
}

bool SlidingMap::add(const Hash &h, int BIT, int FULL) // bit: var with only one bit set
{
	auto it = storage.lower_bound(h);
	
	bool inserted = false;
	if (it != storage.end() && it->first == h) { 
		if (it->second & BIT) return false; 
		it->second |= BIT;
	} else {
		it = storage.insert({h, BIT}).first;
		inserted = true;
	}

	assert(it != storage.end());
	assert(boundary != storage.end() || !query_size);
	// dprn("> Add:{:5} > boundary {} is at position {}", storage.size(), print_it(boundary), distance(storage.begin(), boundary));
	if (query_size && it->first < boundary->first) {
		intersection += (it->second == FULL);
		if (inserted) {
			intersection -= (boundary->second == FULL);
			assert(boundary != storage.begin()); // |S| >= 1!
			boundary--;
			// dprn("> Add:{:5} > boundary {} is at position {} after--", storage.size(), print_it(boundary), distance(storage.begin(), boundary));
		}
	}
	return true;
}

bool SlidingMap::remove(const Hash &h, int BIT, int FULL) // bit: var with only one bit set
{	
	auto it = storage.lower_bound(h); 
	if (it->first != h || !(it->second & BIT)) 
		return false;

	assert(boundary != storage.end() || !query_size);
	// dprn("> Remove:{:5} > boundary {} is at position {}", storage.size(), print_it(boundary), distance(storage.begin(), boundary));
	if (query_size && it->first <= boundary->first) {
		intersection -= (it->second == FULL);
		if (it->second == BIT) {
			boundary++;
			// dprn("> Remove:{:5} > boundary {} is at position {} after++", storage.size(), print_it(boundary), distance(storage.begin(), boundary));
			// assert(boundary != storage.end());
			if (boundary != storage.end())
				intersection += (boundary->second == FULL);
		}
	}

	assert(it != storage.end());
	if (it->second == BIT) {
		assert(it != boundary);
		storage.erase(it);
	} else {
		it->second &= ~BIT;
	}
	return true;
}

void SlidingMap::add_to_query(const Hash &h) 
{	
	if (!add(h, 1)) 
		return;

	limit = relaxed_jaccard_estimate(++query_size, kmer_size) ;
	assert(boundary != storage.end() || (query_size == 1));
	// dprn("> AddQuery:{:5} > boundary {} is at position {}", storage.size(), print_it(boundary), distance(storage.begin(), boundary));
	if (boundary == storage.end()) {
		boundary = storage.begin();
	} else {
		assert(boundary != storage.end());
		boundary++;
	}
	// dprn("> AddQuery:{:5} > boundary {} is at position {} after++", storage.size(), print_it(boundary), distance(storage.begin(), boundary));
	
	assert(boundary != storage.end());
	intersection += (boundary->second == 3);
}

void SlidingMap::remove_from_query(const Hash &h) 
{
	if (!remove(h, 1)) 
		return;
	
	limit = relaxed_jaccard_estimate(--query_size, kmer_size) ;
	assert(boundary != storage.begin() || !query_size);
	// dprn("> RemoveQuery:{:5} > boundary {} is at position {}", storage.size(), print_it(boundary), distance(storage.begin(), boundary));
	if (boundary != storage.end())
		intersection -= (boundary->second == 3);
	if (boundary == storage.begin()) {
		boundary = storage.end();
	} else { 
		boundary--;
	}
	// dprn("> RemoveQuery:{:5} > boundary {} is at position {} after--", storage.size(), print_it(boundary), distance(storage.begin(), boundary));
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
