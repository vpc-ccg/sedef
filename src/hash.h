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

struct Hash {
	enum Status {
		HAS_UPPERCASE, ALL_LOWERCASE, HAS_N
	};
	uint32_t hash; 
	Status status; // 2 if N, 1 if lowercase, 0 otherwise
};

struct Minimizer {
	Hash hash;
	int loc;
};

namespace std {
	template <> struct hash<Hash> {
		size_t operator()(const Hash &h) const
		{
			return std::hash<uint32_t>()(h.hash) ^ std::hash<char>()((char)h.status);
		}
	};
}

/******************************************************************************/

struct Sequence {
	string name;
	string seq;
	const bool is_rc;

	Sequence(const string &name, const string &seq, bool is_rc = false);
};

struct Index {
	const int kmer_size;
	const int window_size;
	unsigned int threshold;

	shared_ptr<Sequence> seq;

	// (hash, loci), sorted by loci
	vector<Minimizer> minimizers;
	// hash -> list of locations
	unordered_map<Hash, list<int>> index;

public:
	Index(shared_ptr<Sequence> seq, int kmer_size = KMER_SIZE, int window_size = WINDOW_SIZE);

	// Find first minimizer at loci p
	int find_minimizers(int p) const;
};

#include "extern/ostream.h"

ostream& operator<<(ostream& os, const Hash& dt);

bool operator<(const Hash &x, const Hash &y);
bool operator<=(const Hash &x, const Hash &y);
bool operator==(const Hash &x, const Hash &y);
bool operator!=(const Hash &x, const Hash &y);
bool operator==(const Minimizer &x, const Minimizer &y);
