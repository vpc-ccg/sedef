/// 786

#pragma once

#include <map>
#include <set>
#include <vector>
#include <string>
#include <queue>
#include "hash.h"

// Keeps (hash, insertion_time) inside
struct SlidingMap: public map<pair<hash_t, int64_t>, char> {

public:
    typename std::map<pair<hash_t, int64_t>, char>::iterator boundary;
    int query_size;
    int intersection;
    int64_t time;

public:
    // sets: query, ref
    SlidingMap(): query_size(0), time(0), intersection(0) { 
        boundary = this->end();
    } 

    void rewind()
    {
        boundary = this->begin();
        std::advance(boundary, query_size - 1);
        assert(boundary != this->end());
        intersection = (boundary->second == 3);
        for (auto it = this->begin(); it != boundary; it++)
            intersection += (it->second == 3);
    }

    int jaccard() 
    {
        if (intersection >= tau() * query_size) {
            return intersection;
        } else {
            return -intersection;
        }
    }

    // &1 = query; &2 = reference

    void add_to_query(const hash_t &h) 
    {
        query_size++;
        if (boundary != this->end() && next(boundary) != this->end()) {
            boundary++;
            intersection += (boundary->second == 3);
        }
        if (h.first) return;

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

    void add_to_reference(const hash_t &h)
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

        if (boundary != this->end() && it->first <= boundary->first) {
            intersection += (it->second == 3);
            if (inserted) {
                intersection -= (boundary->second == 3);
                boundary--;
            }
        }
    }

    void remove_from_reference(const hash_t &h)
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
};


struct SlidingMap2: public multimap<hash_t, bool> {
    typename std::multimap<hash_t, bool>::iterator boundary;
    typename std::multimap<hash_t, bool>::iterator tmp;

public:
    multiset<hash_t> query;
    multiset<hash_t> ref;
    int query_n;

public:
    // sets: query, ref
    SlidingMap2(): query_n(0) { 
        boundary = this->end();
    } 

    void rewind()
    {}

    int jaccard() 
    {
        auto i = query.begin();
        auto j = ref.begin();

        int s = query.size() + query_n;
        int total = 0;
        int jaccard = 0;
        while (i != query.end() && j != ref.end() && total < s) {
            if (*i == *j) { // handle N-s?
                jaccard++; 
                i++; j++;
            } else if (*i < *j) {
                i++;
            } else {
                j++;
            }
            total++;
        }

        // eprn("calc: jaccard={} s={} need={}", jaccard, s, tau()*s);

        if (jaccard >= tau() * s) {
            return jaccard;
        } else {
            return -jaccard;
        }

        // return jaccard() >= tau() * (query.size() + query_n.size());
    }

    void add_to_query(const hash_t &h) 
    {
        if (h.first) {
            query_n++;
            return;
        }
        query.insert(h);
    }

    void add_to_reference(const hash_t &h)
    {
        if (h.first) return;
        ref.insert(h);
    }

    void remove_from_reference(const hash_t &h)
    {
        if (h.first) return;
        auto it = ref.find(h);
        if (it != ref.end()) 
            ref.erase(it);
    }
};

struct Hit {
    int p, q; // query range
    int i, j; // reference range
    double id; // identity score
    /*char break_criteria; /* criteria which prevented further extends 
        0: overlap
        1: trailing Ns
        2: extend right failed */
    int jaccard; // coordinates of seed matches
};

vector<Hit> search (int query_start, 
                    const Hash &ref_hash, 
                    const Hash &query_hash, 
                    bool allow_overlaps=false);

