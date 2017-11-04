/// 786

#include <bits/stdc++.h>
#include <fmt/format.h>
using namespace std;

#define prn(f, ...)    fmt::print(f "\n", __VA_ARGS__)
#define prnn(...)      fmt::print(__VA_ARGS__)

#define eprn(f, ...)   fmt::print(stderr, f "\n", __VA_ARGS__)
#define eprnn(...)     fmt::print(stderr, __VA_ARGS__)

static char dna_lookup[128] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,
    3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,
    0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0
};
inline char qdna(char c) 
{
    return dna_lookup[c];
}


#include <boost/icl/interval_map.hpp>
#include <boost/icl/interval_set.hpp>

typedef boost::icl::discrete_interval<int> INTERVAL;
typedef boost::icl::interval_map<int, boost::icl::interval_map<int, set<int> > >::const_iterator TREE_iterator;
boost::icl::interval_map<int, 
    boost::icl::interval_map<int, set<pair<INTERVAL, INTERVAL> > > > TREE;

int add_pair(INTERVAL a, INTERVAL b) {
    set<pair<INTERVAL, INTERVAL> > s;
    s.insert(make_pair(a, b));
    boost::icl::interval_map<int, set<pair<INTERVAL, INTERVAL> > > iss;
    iss += make_pair(b, s);
    TREE += make_pair(a, iss);
}

int main(void) {
    add_pair(INTERVAL(0, 10), INTERVAL(100, 110));
    add_pair(INTERVAL(5, 15), INTERVAL(105, 115));
    TREE -= INTERVAL(0, 10);

    for (auto &it: TREE) {
        prnn("[{}, {}): {{", it.first.lower(), it.first.upper());
        for (auto &itt: it.second) {
            prnn("[{}, {})->\"", itt.first.lower(), itt.first.upper());
            for (auto &ist: itt.second) prnn("{}-{}:{}-{} ", ist.first.lower(), ist.first.upper(), 
                ist.second.lower(), ist.second.upper());
            prnn("\" ");
        }
        prnn("\n");
    }

    return 0;
}




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
            auto cur = q.front(); q.pop();
            for (int i = 0; i < 4; i++) {
                auto t = cur->child[i];
                if (t) {
                    auto f = cur->fail;
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
            auto cur = q.front(); q.pop();
            for (int i = 0; i < 4; i++) {
                if (cur->child[i]) 
                    q.push(cur->child[i]);
                auto c = cur;
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
    AHOAutomata (): nodes_count(1)
    {

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

                // cnts[ln]++;
                insert_pattern(trie, 0, patterns.size() - 1);
            }
        }
        // PATTERN_PROB = 0;
        // for (map<int, int>::iterator i=cnts.begin();i!=cnts.end();i++)
        //     PATTERN_PROB += (i->second/double(patterns.size())) * pow(1-MAX_EDIT_ERROR, i->first);
        initialize_automata();
    }

    ~AHOAutomata () 
    {
        return;
        queue<AHOTrie*> q;
        q.push(trie);
        while (!q.empty()) {
            auto cur = q.front(); q.pop();
            for (int i = 0; i < 4; i++) {
                if (cur->child[i])
                    q.push(cur->child[i]);
            }
            delete cur;
        }
    }

public:
    void search (const char *text, int len, map<int, int> &hits, int flag) 
    {
        set<int> used;
        AHOTrie *cur = trie, *largest = 0;
        for (int i = 0; i < len; i++) {
            cur = cur->child[qdna(text[i])];
            if (!cur->next_to_output) continue;
            auto x = cur->next_to_output;
            hits[x->output] |= flag;
        }
    }

    string pattern (int i) {
        return patterns[i];
    }
};

int main2(int argc, char **argv)
{
    std::ios::sync_with_stdio(false);
    setlocale(LC_NUMERIC, "en_US.UTF-8");

    AHOAutomata aho;

    string r, q, wr, wq;
    vector<int> qgram_p(1<<(2*8),0), qgram_r(1<<(2*8),0);
    while (cin >> wr >> r >> wq >> q) {
        prnn("{}\t{}\t", wr, wq);
        int common = 0;
        map<int, int> hits;
        aho.search(q.c_str(), q.size(), hits, 1);
        aho.search(r.c_str(), r.size(), hits, 2);
        for (auto &it: hits)
            if (it.second == 3) common++;
        prnn("{}\t", common);

        int maxlen = max(q.size(), r.size());
        int QG = 5;
        uint64_t QSZ = (1 << (2 * QG)); 
        uint64_t MASK = QSZ - 1;
        
        for (uint64_t qi = 0, qgram = 0; qi < q.size(); qi++) {
            qgram = ((qgram << 2) | qdna(q[qi])) & MASK;
            if (qi >= QG - 1) qgram_p[qgram] += 1;
        }
        for (uint64_t qi = 0, qgram = 0; qi < r.size(); qi++) {
            qgram = ((qgram << 2) | qdna(r[qi])) & MASK;
            if (qi >= QG - 1) qgram_r[qgram] += 1;
        }
        int dist = 0;
        for (uint64_t qi = 0; qi < QSZ; qi++)  {
            dist += min(qgram_p[qi], qgram_r[qi]);
            qgram_p[qi] = qgram_r[qi] = 0;
        }
        prnn("{}\t", dist);

        auto patterns = vector<tuple<string, int, int>> {
            make_tuple("#####....#", 6, 4), 
            make_tuple("#####...##", 7, 3),
        };
        for (int pi = 0; pi < patterns.size(); pi++) {
            uint64_t QSZ = (1 << (2 * QG)); 
            uint64_t MASK = QSZ - 1;

            auto &pattern = get<0>(patterns[pi]);
            for (uint64_t qi = 0, qgram = 0; qi < q.size() - pattern.size(); qi++) {
                for (int i = 0; i < pattern.size(); i++) if (pattern[i] != '.')
                    qgram = ((qgram << 2) | qdna(q[qi + i])) & MASK;
                qgram_p[qgram] += 1;
            }
            for (uint64_t qi = 0, qgram = 0; qi < r.size() - pattern.size(); qi++) {
                for (int i = 0; i < pattern.size(); i++) if (pattern[i] != '.')
                    qgram = ((qgram << 2) | qdna(r[qi + i])) & MASK;
                qgram_r[qgram] += 1;
            }
            int qdist = 0;
            for (uint64_t qi = 0; qi < QSZ; qi++)  {
                qdist += min(qgram_p[qi], qgram_r[qi]);
                qgram_p[qi] = qgram_r[qi] = 0;
            }
            prnn("{}\t", qdist);
        }
        prnn("\n");
    }

    return 0;
}