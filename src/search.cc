/// 786
/// Based on http://www.biorxiv.org/content/biorxiv/early/2017/03/24/103812.full.pdf

#include <list>
#include <queue>
#include <vector>
#include <string>
#include <cmath>
#include <queue>
#include "common.h"
#include "search.h"
using namespace std;

#include <boost/icl/interval_map.hpp>
#include <boost/icl/interval_set.hpp>


typedef boost::icl::discrete_interval<int> INTERVAL;
typedef boost::icl::interval_map<int, set<pair<INTERVAL, INTERVAL> > > SUBTREE_t;
typedef boost::icl::interval_map<int, SUBTREE_t> TREE_t;
TREE_t TREE;


const int MAX_MATCH = 1 * 500 * 1000; // 0.5 MB at most!
const int RIGHT_ALLOWANCE = 750; /// TODO more mathy formulation

int QGRAM_NORMAL_FAILED = 0;
int QGRAM_SPACED_FAILED = 0;
int CORE_FAILED = 0;


inline bool in_map_n(const map<hash_t, bool> &m, const hash_t &k) { 
    return !k.first && m.find(k) != m.end();
}

inline int min_qgram(int l, int q) {
    return l * (1 - MAX_GAP_ERROR - q * MAX_EDIT_ERROR) - (GAP_FREQUENCY * l + 1) * (q - 1);
}

inline bool hits_cmp(const Hit &a, const Hit &b) {
    if (a.q != b.q) return a.q > b.q; // .p are all equal
    return make_pair(a.i, a.j) < make_pair(a.i, a.j);   // larger intervals at the top! then by other ints
}

inline bool check_overlap (TREE_t::const_iterator &pf, int pf_pos, int pfp_pos) {
    SUBTREE_t::const_iterator pfp;
    if (pf != TREE.end() && (pfp = pf->second.find(pfp_pos)) != pf->second.end()) { // check overlap!
        for (set<pair<INTERVAL, INTERVAL> >::iterator it = pfp->second.begin(); it != pfp->second.end(); it++) {
            int sA = it->first.lower(),  eA = it->first.upper();
            int sB = it->second.lower(), eB = it->second.upper();
            assert(eA - sA == eB - sB);
            if (eA - sA > MIN_READ_SIZE * 1.5 && eA - pf_pos >= RIGHT_ALLOWANCE && eB - pfp_pos >= RIGHT_ALLOWANCE)
                return false;
        }
    }
    return true;
}

inline void add_to_tree(INTERVAL a, INTERVAL b) {
    set<pair<INTERVAL, INTERVAL> > s;
    s.insert(make_pair(a, b));
    boost::icl::interval_map<int, set<pair<INTERVAL, INTERVAL> > > iss;
    iss += make_pair(b, s);
    TREE += make_pair(a, iss);
}


AHOAutomata *aho = NULL;
vector<int> qgram_p, qgram_r;

pair<int, int> filter(const string &q, int q_pos, int q_len, const string &r, int r_pos, int r_len) 
{
    int maxlen = max(q_len, r_len);
    int QG = 5;
    uint64_t QSZ = (1 << (2 * QG)); 
    uint64_t MASK = QSZ - 1;

    int minqg = min_qgram(maxlen, QG);
    assert(minqg >= 10);
    
    for (uint64_t qi = q_pos, qgram = 0; qi < q_pos + q_len; qi++) {
        qgram = ((qgram << 2) | qdna(q[qi])) & MASK;
        if (qi - q_pos >= QG - 1) qgram_p[qgram] += 1;
    }
    for (uint64_t qi = r_pos, qgram = 0; qi < r_pos + r_len; qi++) {
        qgram = ((qgram << 2) | qdna(r[qi])) & MASK;
        if (qi - r_pos >= QG - 1) qgram_r[qgram] += 1;
    }
    int dist = 0;
    for (uint64_t qi = 0; qi < QSZ; qi++)  {
        dist += min(qgram_p[qi], qgram_r[qi]);
        qgram_p[qi] = qgram_r[qi] = 0;
    }

    if (dist < minqg) {
        QGRAM_NORMAL_FAILED++;
        return make_pair(-1, -1);
    }
   
    // w=50; k=5; s=11; k=8 -->7
    // auto patterns = vector<tuple<string, int, int>> {
    //     make_tuple("#####....#", 6, 4), 
    //     make_tuple("#####...##", 7, 3),
    // };
    // for (int pi = 0; pi < patterns.size(); pi++) {
    //     int minqg = min_qgram(maxlen, QG);
    //     uint64_t QSZ = (1 << (2 * QG)); 
    //     uint64_t MASK = QSZ - 1;

    //     auto &pattern = get<0>(patterns[pi]);
    //     for (uint64_t qi = q_pos, qgram = 0; qi < q_pos + q_len - pattern.size(); qi++) {
    //         for (int i = 0; i < pattern.size(); i++) if (pattern[i] != '.')
    //             qgram = ((qgram << 2) | qdna(q[qi + i])) & MASK;
    //         qgram_p[qgram] += 1;
    //     }
    //     for (uint64_t qi = r_pos, qgram = 0; qi < r_pos + r_len - pattern.size(); qi++) {
    //         for (int i = 0; i < pattern.size(); i++) if (pattern[i] != '.')
    //             qgram = ((qgram << 2) | qdna(r[qi + i])) & MASK;
    //         qgram_r[qgram] += 1;
    //     }
    //     int qdist = 0;
    //     for (uint64_t qi = 0; qi < QSZ; qi++)  {
    //         qdist += min(qgram_p[qi], qgram_r[qi]);
    //         qgram_p[qi] = qgram_r[qi] = 0;
    //     }
    //     if (qdist < minqg - get<2>(patterns[pi])) { // this should be worked upon a little bit
    //         QGRAM_SPACED_FAILED++;
    //         return make_pair(-1, -1);
    //     }
    // }

    map<int, int> hits;
    assert(q_len == r_len);
    int common = 0;
    aho->search(q.c_str() + q_pos, q_len, hits, 1);
    aho->search(r.c_str() + r_pos, r_len, hits, 2);
    for (map<int, int>::iterator it = hits.begin(); it != hits.end(); it++)
        if (it->second == 3) common++;

    double boundary = (3.0/4) * (q_len / 50.0);
    if (common < boundary) {
        CORE_FAILED++;
        return make_pair(-1, -1);
    }

    return make_pair(dist, common);
}

void refine(int p, int idx_p, int idx_q, // query start and W(query) range
            int x, int y, // window to check for the initial MIN_READ_SIZE match
            int initial_match_length,
            //map<hash_t, char> 
            map<hash_t, bool> L0, // W(query)
            const Hash &ref_hash, 
            const Hash &query_hash,
            vector<Hit> &hits, // Output
            bool allow_overlaps=false)
{
    double tau = ::tau();

    // eprn(">> eval {} vs {}", p, x);

    TREE_t::const_iterator pf = TREE.find(p);
    if (!check_overlap(pf, p, x)) {
        return;
    }
    
    int i = x, 
        j = x + initial_match_length;
    // L0 is |W(query)| and L0[h] is {location of h in the query, false}
    // L  is |W(reference) + W(query)| 
    // L[h] is {location of h in the reference, true} iff h is in |W(reference) & W(query)|
    SlidingMap<hash_t, bool> L(L0); // matches: ref -> query!
    int idx_i = ref_hash.find_minimizers(i);
    if (idx_i == ref_hash.minimizers.size()) {
        return;
    }
    int idx_j = idx_i;
    while (idx_j < ref_hash.minimizers.size() && ref_hash.minimizers[idx_j].second < j) {
        // do not set 1 to N hashes
        L.store[ref_hash.minimizers[idx_j].first] = in_map_n(L0, ref_hash.minimizers[idx_j].first);
        idx_j++;
    }

    L.boundary = L.store.begin();
    std::advance(L.boundary, L0.size() - 1);
    
    int jaccard = L.boundary->second;
    for (map<hash_t, bool>::iterator it = L.store.begin(); it != L.boundary; it++)
        jaccard += it->second;

    // Extend it up to MIN_READ_SIZE (initial match)
    int seed_i = -1, seed_jaccard = 0;
    if (jaccard >= tau * L0.size())
        seed_i = i, seed_jaccard = jaccard;
    // prn("[{} {}] ja={} se_i={} se_ja={}", x, y, jaccard, seed_i, seed_jaccard);
    while (i < y) {
        if (idx_i < ref_hash.minimizers.size() && ref_hash.minimizers[idx_i].second < i + 1) {
            jaccard += L.remove(ref_hash.minimizers[idx_i].first);
            idx_i++;        
        }
        if (idx_j < ref_hash.minimizers.size() && ref_hash.minimizers[idx_j].second < j + 1) {
            jaccard += L.add(ref_hash.minimizers[idx_j].first, in_map_n(L0, ref_hash.minimizers[idx_j].first));
            idx_j++;
        }
        if (jaccard >= tau * L0.size() && jaccard >= seed_jaccard) {
            seed_i = i, seed_jaccard = jaccard;
            break;
        }
        i++, j++;
    }
    if (seed_i == -1) {
        return;
    }

    // Now keep extending it as much as we can
    // L has all minimizers now
    // we are mapping query[p, q] to ref[i, j]
    int q = p + initial_match_length;
    // idx_i/idx_j points to the first/last minimizer of REF_i,j
    // idx_p/idx_q points to the first/last minimizer of QUERY_p,q
    double prev_jaccard = double(jaccard) / L0.size();
    double init_jaccard = prev_jaccard;
    int prev_jaccard_p = jaccard;
    int init_jaccard_p = jaccard;
    
    const int RECOVER_BP = 250; // We allow 250bp extra extend just in case!
    int recover = 0;

    int best_q = q;
    int best_j = j;
    char reason = 0;
    while (q < query_hash.seq.size() && j < ref_hash.seq.size()) {    
        // - disallow overlaps or too long matches
        int max_match = MAX_MATCH;
        if (!allow_overlaps) 
            max_match = min(max_match, int((1.0 / MAX_GAP_ERROR + .5) * abs(p - i)));
        if (max(q - p, j - i) > max_match) {
            reason = 0;
            break;
        }

        // - try extending to the right
        if (query_hash.minimizers[idx_q].second < q + 1) { // extend query
            hash_t h = query_hash.minimizers[idx_q].first;
            if (!in_map(L0, h)) {
                L0[h] = 0; // {query_hash.minimizers[idx_q].second, false};
                // check is h in L now
                jaccard += L.add(h, in_map_n(L.store, h));
                L.boundary++; // |W(A)| increases
                jaccard += L.boundary->second;
            }
            idx_q++;
        }
        q++;
        if (idx_j < ref_hash.minimizers.size() && ref_hash.minimizers[idx_j].second < j + 1) { // extend reference
            const minimizer_t &h = ref_hash.minimizers[idx_j];
            jaccard += L.add(h.first, in_map_n(L0, h.first));
            idx_j++;
        } 
        j++;
        if (jaccard >= tau * L0.size()) {
            prev_jaccard = double(jaccard) / L0.size();
            prev_jaccard_p = jaccard;
            best_j = j;
            best_q = q;
            recover = 0;
            continue;
        } 

        if (recover < RECOVER_BP) {
            recover++;
            continue;
        }
        reason = 2;
        break;
    }

    pair<int, int> edist = filter(query_hash.seq, p, best_q - p + 1, ref_hash.seq, i, best_j - i + 1);
    if (edist.first < 0)
        return;

    hits.push_back(Hit(p, best_q, i, best_j, j2md(prev_jaccard), j2md(init_jaccard), reason, make_pair(init_jaccard_p, prev_jaccard_p), edist));
    add_to_tree(INTERVAL(p, best_q), INTERVAL(i, best_j));
}

vector<Hit> search (int query_start, 
                    const Hash &ref_hash, 
                    const Hash &query_hash, 
                    bool allow_overlaps) 
{

    if (!qgram_p.size()) {
        qgram_p = vector<int>(1<<(2*8),0);
        qgram_r = vector<int>(1<<(2*8),0);
        aho = new AHOAutomata();
    }

    map<hash_t, bool> L0;
    int st = query_hash.find_minimizers(query_start);
    if (st == query_hash.minimizers.size())
        return vector<Hit>();
    int mi = st;
    //prn("{}\n", mi);
    
    TREE_t::const_iterator pf = TREE.find(query_start);

    // Iterate through all unique hashes in query[query_start: query_start + MIN_READ_SIZE]
    // `candidates` is a list of positions in the reference which match the query hashes
    vector<int> candidates;
    for (; mi < query_hash.minimizers.size() && query_hash.minimizers[mi].second - query_start <= MIN_READ_SIZE; mi++) { 
        const minimizer_t &m = query_hash.minimizers[mi];
        if (mi != st && m.first == query_hash.minimizers[mi - 1].first) 
            continue;
        
        L0[m.first] = false;
        // If it is N hash, ignore it
        if (m.first.first)
            continue;
        // Is this hash in the reference? Is this hash high-call hash as well?
        map<hash_t, list<int>, MapCompare>::const_iterator ptr = ref_hash.index.find(m.first);
        if (ptr == ref_hash.index.end() || ptr->second.size() >= ref_hash.threshold) 
            continue;
        // Iterate through positions in the reference with this hash
        for (list<int>::const_iterator PI = ptr->second.begin(); PI != ptr->second.end(); PI++) {
            int pos = *(PI);
            // Make sure to have at least 1 kb spacing if reference = query
            if (!allow_overlaps && pos < query_start + 2 * MIN_READ_SIZE)
                continue;
            if (!check_overlap(pf, query_start, pos)) 
                continue;
            candidates.push_back(pos);
        }
    }
    sort(candidates.begin(), candidates.end());
    
    if (L0.size() < 1) { // TODO min sketch size limit for Ns
        return vector<Hit>();
    }
    int M = ceil(L0.size() * tau());
    // eprn("{} candidates to evaluate (M={}, L0={}, tau={})", candidates.size(), M, L0.size(), tau());

    // Find all locations in the `candidates` so that
    // MIN_READ_SIZE read covers at least M (= s * tau) hashes
    vector<pair<int, int> > T; T.push_back(make_pair(-1,-1));
    for (int i = 0; i <= (int)candidates.size() - M; i++) {
        int j = i + (M - 1);

        if (candidates[j] - candidates[i] < MIN_READ_SIZE) {
            int x = max(0, candidates[j] - MIN_READ_SIZE + 1), 
                y = candidates[i] + 1;
            if (x >= T.back().second)
                T.push_back(make_pair(x, y));
            T.back().second = y;
        }
    }

    vector<Hit> hits, hits_real;
    for (int t = 1; t < T.size(); t++) {
        refine(
            query_start, st, mi, 
            T[t].first, T[t].second, 
            MIN_READ_SIZE,
            L0, ref_hash, query_hash, hits, allow_overlaps
        );
    }
    sort(hits.begin(), hits.end(), hits_cmp);
    for (int HI = 0; HI < hits.size(); HI++) { // brute force, but maybe easier than doing full interval thing
        const Hit &h = hits[HI];
        bool add = true;
        for (int j = hits_real.size() - 1; j >= 0; j--) {
            const Hit &ph = hits_real[j];
            if (h.i >= ph.i && h.j <= ph.j) { // if full match
                add = false;
                break;
            }
        }
        if (add) hits_real.push_back(h);
    }
    reverse(hits_real.begin(), hits_real.end()); // smaller to larger
    TREE -= INTERVAL(0, query_start);

    return hits_real;
}

