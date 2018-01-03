/// 786

/******************************************************************************/

#pragma once

/******************************************************************************/

#include <list>
#include <iostream>
#include <unordered_map>
#include <vector>

using namespace std;

/******************************************************************************/

typedef pair<char, uint32_t> hash_t; // 2 if N, 1 if lowercase, 0 otherwise
typedef pair<hash_t, int> minimizer_t;

/******************************************************************************/

struct Hash {
	const int kmer_size;
	const int window_size;

	unsigned int threshold;
	string seq;

	// (hash, loci), sorted by loci
	vector<minimizer_t> minimizers;
	// hash -> list of locations
	unordered_map<hash_t, list<int>, pair_hash> index;

public:
	Hash(const string &s, int kmer_size = KMER_SIZE, int window_size = WINDOW_SIZE);

	// Find first minimizer at loci p
	int find_minimizers(int p) const;
};


#include "extern/ostream.h"
ostream& operator<<(ostream& os, const hash_t& dt);
