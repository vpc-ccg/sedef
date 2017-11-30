/// 786

#pragma once

#include <map>
#include <set>
#include <vector>
#include <string>
#include <queue>
#include "hash.h"

template<typename K>
struct SlidingMap: public multimap<K, bool> {
    typename std::multimap<K, bool>::iterator boundary;
    typename std::multimap<K, bool>::iterator tmp;

    inline auto match_winnow_query(multimap<K, bool> *m, const K &k) 
    { 
        if (k.first) return make_pair(false, m->end());
        for (auto it = m->lower_bound(k); it != m->end() && it->first == k; it++) {
            if (!it->second) return make_pair(true, it);
        }
        return make_pair(false, m->end());
    }

public:
    multimap<K, bool> query;
    int jaccard;

public:
    // sets: query, ref
    SlidingMap(): jaccard(0) {} 

    int add(const K &key, bool val) 
    {
        int diff = 0;
        if (key < boundary->first) {
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
        std::advance(boundary, query.size() - 1);
        
        jaccard = boundary->second;
        for (auto it = this->begin(); it != boundary; it++)
            jaccard += it->second;
    }

    bool valid_jaccard() 
    {
        return jaccard >= tau() * query.size();
    }

    int add_to_query(const K &h) 
    {
        auto it = query.insert({h, false});

        auto is_in = match_winnow_query(this, h);
        if (is_in.first) {
            it->second = is_in.second->second = true; // mark hash in both winnows as matched!
        }
        jaccard += this->add(h, is_in.first);
        boundary++; // |W(A)| increases
        jaccard += boundary->second;
    }

    void add_to_reference(const K &h)
    {
        auto is_in = match_winnow_query(&query, h);
        auto it = this->insert({h, is_in.first});
        if (is_in.first) {
            is_in.second->second = true; // it is marked as NOT TOUCHED [fix: pointer maybe?]
        }
    }

    void remove_from_reference(const K &h)
    {
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
    double id, init_id; // identity score
    char break_criteria; /* criteria which prevented further extends 
        0: overlap
        1: trailing Ns
        2: extend right failed */
    pair<int, int> jaccard; // coordinates of seed matches
    pair<int, int> edist;

    Hit(int p, int q, int i, int j, double id, double init_id, char break_criteria, pair<int, int> jaccard, pair<int, int> edist):
        p(p), q(q), i(i), j(j), id(id), init_id(init_id), break_criteria(break_criteria), jaccard(jaccard), edist(edist) {}
};

vector<Hit> search (int query_start, 
                    const Hash &ref_hash, 
                    const Hash &query_hash, 
                    bool allow_overlaps=false);

extern char _binary_patterns_bin_start;
extern char _binary_patterns_bin_end;
extern char _binary_patterns_bin_size;

class AHOAutomata {
public:
    struct AHOTrie {
        AHOTrie *child[4];        
        AHOTrie *fail;            
        AHOTrie *next_to_output;  
        int level;                
        int id;                   
        uint64_t bin_size;        
        int32_t output;           
        list<void*> bin;

    public:
        AHOTrie ():
            fail(0), output(-1), bin_size(0), next_to_output(0), level(0), id(0)
        {
            memset(child, 0, 4 * sizeof(AHOTrie*));
        }
    };

private:
    vector<string> patterns;
    AHOTrie *trie;
    int nodes_count;

private:
    void insert_pattern (AHOTrie *n, int level, int id) 
    {
        if (level < patterns[id].size()) {
            char cx = qdna(patterns[id][level]);
            if (n->child[cx] == 0) {
                n->child[cx] = new AHOTrie();
                n->child[cx]->level = level + 1;
                nodes_count++;
            }
            insert_pattern(n->child[cx], level + 1, id); 
        }
        else n->output = id;
    }

    void initialize_automata () 
    {
        queue<AHOTrie*> q;
        trie->fail = trie;

        for (int i = 0; i < 4; i++) 
            if (trie->child[i]) {
                trie->child[i]->fail = trie;
                q.push(trie->child[i]);
            }

        int traversed = 0;
        while (!q.empty()) {
            AHOTrie *cur = q.front(); q.pop();
            for (int i = 0; i < 4; i++) {
                AHOTrie *t = cur->child[i];
                if (t) {
                    AHOTrie *f = cur->fail;
                    while (f != trie && !f->child[i])
                        f = f->fail;
                    t->fail = (f->child[i] ? f->child[i] : trie);
                    q.push(t);
                    t->next_to_output = (t->fail->output >= 0) ? t->fail : t->fail->next_to_output;
                }
            }
            cur->id = ++traversed;
        }
        q.push(trie);
        while (!q.empty()) {
            AHOTrie *cur = q.front(); q.pop();
            for (int i = 0; i < 4; i++) {
                if (cur->child[i]) 
                    q.push(cur->child[i]);
                AHOTrie *c = cur;
                while (c != trie && !c->child[i]) 
                    c = c->fail;
                cur->child[i] = (c->child[i] ? c->child[i] : trie);
                c = cur;
                while (c && c->output == -1)
                    c = c->next_to_output;
                cur->next_to_output = c;
            }
        }
    }

public:
    double PATTERN_PROB;
    AHOAutomata (): nodes_count(1)
    {
        map<int, int> cnts;

        trie = new AHOTrie();
        patterns.reserve(5000000);

        char *data = &_binary_patterns_bin_start;
        int pos = 0;
        int64_t size = (int64_t)(&_binary_patterns_bin_size);
        while (1) {
            int16_t ln;
            if (pos == size)
                break;
            memcpy(&ln, data + pos, sizeof(int16_t));
            pos += sizeof(int16_t);
            
            int32_t cnt;
            memcpy (&cnt, data + pos, sizeof(int32_t));
            pos += sizeof(int32_t);
            
            int sz = ln / 4 + (ln % 4 != 0);
            int64_t x;
            for(int i = 0; i < cnt; i++) {
                memcpy (&x, data + pos, sz);
                assert(sz <= 8);
                pos += sz;
                
                patterns.push_back("");
                for(int j = ln - 1; j >= 0; j--)
                    patterns.back() += "ACGT"[(x >> (2 * j)) & 3];

                if (patterns.back() == "AAAAAAAAAAAAAA") continue;
                if (patterns.back() == "CCCCCCCCCCCCCC") continue;
                if (patterns.back() == "GGGGGGGGGGGGGG") continue;
                if (patterns.back() == "TTTTTTTTTTTTTT") continue;

                cnts[ln]++;
                insert_pattern(trie, 0, patterns.size() - 1);
            }
        }
        PATTERN_PROB = 0;
        for (map<int, int>::iterator i = cnts.begin(); i != cnts.end(); i++)
            PATTERN_PROB += (i->second / double(patterns.size())) * pow(1 - MAX_EDIT_ERROR, i->first);
        // sort(patterns.begin(), patterns.end());
        // for (int i = 0; i < patterns.size(); i++)
        //     prn("{}", patterns[i]);
        // exit(0);
        initialize_automata();
    }

public:
    void search (const char *text, int len, map<int, int> &hits, int flag) 
    {
        set<int> used;
        AHOTrie *cur = trie, *largest = 0;
        for (int i = 0; i < len; i++) {
            cur = cur->child[qdna(text[i])];
            if (!cur->next_to_output) continue;
            AHOTrie *x = cur->next_to_output;
            hits[x->output] |= flag;
        }
    }

    string pattern (int i) {
        return patterns[i];
    }
};