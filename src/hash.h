/// 786

#pragma once

#include <list>
#include <vector>
using namespace std;

typedef pair<bool, uint32_t> hash_t; // 1 if does not contain N, 0 otherwise
typedef pair<hash_t, int> minimizer_t;

struct MapCompare {
    bool operator() (const hash_t& lhs, const hash_t& rhs) const {
        if (lhs.second != rhs.second) return lhs.second < rhs.second;
        return lhs.first < rhs.first;
    }
};

struct Hash {
    unsigned int threshold;
    string seq;

    // (hash, loci), sorted by loci
    vector<minimizer_t> minimizers;
    // hash -> list of locations
    map<hash_t, list<int>, MapCompare> index;

public:
    Hash() {};
    Hash(const string &s);
    // Find first minimizer at loci p
    int find_minimizers(int p) const;
};
