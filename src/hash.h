/// 786

#pragma once

#include <list>
#include <vector>
#include <unordered_map>
using namespace std;

typedef pair<bool, uint32_t> hash_t; // 1 if does not contain N, 0 otherwise
typedef pair<hash_t, int> minimizer_t;

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const {
        return hash<T2>{}(p.second);
    }
};

typedef array<int, 4> qgram_t;

struct Hash {
    unsigned int threshold;
    string seq;

    vector<qgram_t> qgram;

    // (hash, loci), sorted by loci
    vector<minimizer_t> minimizers;
    // hash -> list of locations
    unordered_map<hash_t, list<int>, pair_hash> index;

public:
    Hash() {};
    Hash(const string &s);
    // Find first minimizer at loci p
    int find_minimizers(int p) const;
};
