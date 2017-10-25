/// 786

#pragma once
#include <map>
#include "hash.h"

template<typename K, typename V>
struct SlidingMap {
    std::map<K, V> store;
    typename std::map<K, V>::iterator boundary;

public:
    SlidingMap(std::map<K, V> m): store(m) {}
    SlidingMap() {} 
    

    int add(K key, V val) 
    {
        int diff = 0;
        auto present = store.find(key);
        if (key <= boundary->first) {
            if (present == store.end()) {
                diff = val - boundary->second;
                boundary--;
            } else {
                diff = val - present->second;
            }
        }
        store[key] = val;
        return diff;
    }

    int remove(K key) 
    {
        int diff = 0;
        auto present = store.find(key);
        if (present == store.end()) 
            return diff;
        if (key <= boundary->first) {
            boundary++;
            diff = boundary->second - present->second;
        }
        store.erase(present);
        return diff;
    }
};

struct Hit {
    int p, q; // query range
    int i, j; // reference range
    double id, init_id; // identity score
    char break_criteria; /* criteria which prevented further extends 
        0: overlap
        1: trailing Ns
        2: extend right failed */
    pair<int, int> jaccard; // coordinates of seed matches
    int edist;
};

vector<Hit> search (int query_start, 
                    const Hash &ref_hash, 
                    const Hash &query_hash, 
                    bool allow_overlaps=false);