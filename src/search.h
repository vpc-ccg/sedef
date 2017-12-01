/// 786

#pragma once

#include <map>
#include <set>
#include <vector>
#include <string>
#include <queue>
#include "hash.h"

struct SlidingMap: public multimap<hash_t, bool> {
    typename std::multimap<hash_t, bool>::iterator boundary;
    typename std::multimap<hash_t, bool>::iterator tmp;

    inline auto match_winnow_query(multimap<hash_t, bool> *m, const hash_t &k) 
    { 
        if (k.first) return make_pair(false, m->end());
        for (auto it = m->lower_bound(k); it != m->end() && it->first == k; it++) {
            if (!it->second) return make_pair(true, it);
        }
        return make_pair(false, m->end());
    }

public:
    multimap<hash_t, bool> query;
    int query_n;
    int jaccard;

public:
    // sets: query, ref
    SlidingMap(): jaccard(0), query_n(0) { 
        boundary = this->end();
    } 

    int add(const hash_t &key, bool val) 
    {
        int diff = 0;
        if (boundary != this->end() && key < boundary->first) {
            auto present = this->find(key);
            if (present == this->end()) {
                diff = val - boundary->second;
                boundary--;
            } else {
                diff = val - present->second;
            }
        }
        this->insert(make_pair(key, val)); // assume inserted after last
        return diff;
    }

    void rewind()
    {
        boundary = this->begin();
        std::advance(boundary, (query.size() + query_n.size()) - 1);
        
        jaccard = boundary->second;
        for (auto it = this->begin(); it != boundary; it++)
            jaccard += it->second;
    }

    bool valid_jaccard() 
    {
        return jaccard >= tau() * (query.size() + query_n.size());
    }

    void add_to_query(const hash_t &h) 
    {
        if (h.first) {
            query_n++;
            return;
        }

        auto it = query.insert({h, false});

        auto is_in = match_winnow_query(this, h);
        if (is_in.first) {
            it->second = is_in.second->second = true; // mark hash in both winnows as matched!
        }
        jaccard += this->add(h, is_in.first);
        boundary++; // |W(A)| increases
        jaccard += boundary->second;
    }

    void add_to_reference(const hash_t &h)
    {
        if (h.first) return;

        auto is_in = match_winnow_query(&query, h);
        this->insert({h, is_in.first});
        if (is_in.first) {
            is_in.second->second = true; // it is marked as NOT TOUCHED [fix: pointer maybe?]
        }
    }

    void remove_from_reference(const hash_t &h)
    {
        if (h.first) return;

        // LOWER BOUND: remove matching hash "just in case"
        auto itu = this->lower_bound(h); 
        if (itu == this->end() || itu->first != h) 
            return;
        
        bool in_union = false;
        // Find first item which equals <h> and which is in both sets
        for (auto itu2 = itu; itu2 != this->end() && itu2->first == h; itu2++)
            if (itu2->second) {
                itu = itu2; 
                break;
            }

        // itu is "the hash"
        // mark it zero in query set (if it exists)
        for (auto tmpit = query.lower_bound(h); tmpit != query.end() && tmpit->first == h; tmpit++)
            if (tmpit->second) {
                tmpit->second = false;
                break;
            }

        if (h < boundary->first) {
            boundary++;
            jaccard += boundary->second - itu->second;
        } else if (h == boundary->first) {
            for (auto itu2 = itu; itu2 != this->end() && itu2->first == h; itu2++)
                if (itu2 == boundary) {
                    boundary++;
                    jaccard += boundary->second - itu->second;
                    break;
                }
        }
        this->erase(itu);      
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

