/// 786

#pragma once

#include <list>
#include <vector>
#include <unordered_map>
using namespace std;

typedef pair<uint32_t, int> minimizer_t;

struct Hash {
    unsigned int threshold;
    string seq;

    // (hash, loci), sorted by loci
    vector<minimizer_t> minimizers;
    // hash -> list of locations
    unordered_map<uint32_t, list<int>> index;

public:
    Hash(const string &s);
    // Find first minimizer at loci p
    int find_minimizers(int p) const;
};
