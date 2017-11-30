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
typedef boost::icl::interval_map<int, set<pair<INTERVAL, INTERVAL>>> SUBTREE_t;
typedef boost::icl::interval_map<int, SUBTREE_t> TREE_t;
TREE_t TREE;

// Should we limit this at all?!
const int MAX_MATCH = 500 * 1000; // 50kb at most
const int RIGHT_ALLOWANCE = 750; /// TODO more mathy formulation
const int RECOVER_BP = 10; // We allow 250bp extra extend just in case!


int64_t TOTAL_ATTEMPTED = 0;
int64_t JACCARD_FAILED = 0;
int64_t QGRAM_NORMAL_FAILED = 0;
int64_t OTHER_FAILED = 0;
int64_t CORE_FAILED = 0;
int64_t INTERVAL_FAILED = 0;

inline int min_qgram(int l, int q) 
{
    return l * (1 - MAX_GAP_ERROR - q * MAX_EDIT_ERROR) - (GAP_FREQUENCY * l + 1) * (q - 1);
}

inline bool check_overlap (TREE_t::const_iterator &pf, int pf_pos, int pfp_pos) 
{
    SUBTREE_t::const_iterator pfp;
    if (pf != TREE.end() && (pfp = pf->second.find(pfp_pos)) != pf->second.end()) { // check overlap!
        for (auto &it: pfp->second) {
            int sA = it.first.lower(),  eA = it.first.upper();
            int sB = it.second.lower(), eB = it.second.upper();
            assert(eA - sA == eB - sB);
            if (eA - sA > MIN_READ_SIZE * 1.5 
                // && abs((pf_pos - sA) - (pfp_pos - sB)) >= (eA-sA)/4
                && eA - pf_pos  >= RIGHT_ALLOWANCE 
                && eB - pfp_pos >= RIGHT_ALLOWANCE
            ) return false;
        }
    }
    return true;
}

AHOAutomata *aho = NULL;
vector<int> qgram_p, qgram_r;

const int LOW_LIMIT = 500;

auto filter(const string &q, int q_pos, int q_len, const string &r, int r_pos, int r_len) 
{
    int q_up = 0; for (int i = 0; i < q_len; i++) q_up += isupper(q[q_pos + i]); 
    int r_up = 0; for (int i = 0; i < r_len; i++) r_up += isupper(r[r_pos + i]);

    if (q_up < LOW_LIMIT || r_up < LOW_LIMIT) {
        OTHER_FAILED++;
        return make_pair(make_pair(-1, -1), make_pair(-1, -1));
    }
    
    int maxlen = max(q_len, r_len);
    int QG = 5;
    uint64_t QSZ = (1 << (2 * QG)); 
    uint64_t MASK = QSZ - 1;

    int minqg = min_qgram(maxlen, QG);
    assert(minqg >= 10);

    assert(q_pos + q_len < q.size());
    assert(r_pos + r_len < r.size());
    
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
        return make_pair(make_pair(-1, -1), make_pair(-1, -1));
    }

    map<int, int> hits;
    assert(q_len == r_len);
    int common = 0;
    aho->search(q.c_str() + q_pos, q_len, hits, 1);
    aho->search(r.c_str() + r_pos, r_len, hits, 2);
    for (auto &h: hits) if (h.second == 3) common++;

    double boundary = (1.0/2) * (q_len / 50.0);
    if (common < boundary) {
        CORE_FAILED++;
        return make_pair(make_pair(-1, -1), make_pair(-1, -1));
    }

    return make_pair(make_pair(dist, common), make_pair(q_up, r_up));
}

auto parse_hits(vector<Hit> &hits)
{
    vector<Hit> hits_real;
    sort(hits.begin(), hits.end(), [](const Hit &a, const Hit &b) {
        if (a.q != b.q) return a.q > b.q; // .query_start are all equal
        return make_pair(a.i, a.j) < make_pair(a.i, a.j);   // larger intervals at the top! then by other ints
    });
    for (auto &h: hits) { // brute force, but maybe easier than doing full interval thing
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
    return hits_real;
}

void extend(SlidingMap<hash_t> &winnow,
    const Hash &query_hash, int query_start, int query_end,
    int query_winnow_start, int query_winnow_end,
    const Hash &ref_hash, int ref_start, int ref_end,
    int ref_winnow_start, int ref_winnow_end,
    vector<Hit> &hits,
    bool allow_overlaps) 
{
    // eprn(">> extend query:{:n} vs ref:{:n}...{:n}", query_start, min_ref, max_ref);
            
    auto edist = filter(
        query_hash.seq, query_start, query_end - query_start + 1, 
        ref_hash.seq, ref_start, ref_end - ref_start + 1
    );
    if (edist.first.first < 0) return;

    int qlow = query_end - query_start - edist.second.first;
    int rlow = ref_end - ref_start - edist.second.second;

    int best_query_end = query_end;
    int best_ref_end = ref_end;
    double prev_jaccard = double(winnow.jaccard) / winnow.query.size();
    int prev_jaccard_p = winnow.jaccard;
    for (int recover = 0; query_end < query_hash.seq.size() && ref_end < ref_hash.seq.size() && recover < RECOVER_BP; ) {   
        if (qlow > (query_end-query_start)/4 || rlow > (ref_end-ref_start)/4) { // TODO FIX
            break;
        }
        if (islower(query_hash.seq[query_end])) qlow++;
        if (islower(ref_hash.seq[ref_end])) rlow++;

        // - disallow overlaps or too long matches
        int max_match = MAX_MATCH;
        if (!allow_overlaps) {
            max_match = min(max_match, int((1.0 / MAX_GAP_ERROR + .5) * abs(query_start - ref_start)));
        }
        if (max(query_end - query_start, ref_end - ref_start) > max_match) {
            break;
        }

        // N detection
        int n_count = 0;
        while ((query_hash.seq[query_end] == 'N') || (ref_hash.seq[ref_end] == 'N')) {
            n_count++, query_end++, ref_end++;
            if (n_count > 100) goto E;
        }

        // - try extending to the right
        if (query_hash.minimizers[query_winnow_end].second < query_end + 1) { // extend query
            winnow.add_to_query(query_hash.minimizers[query_winnow_end].first);
            query_winnow_end++;
        }
        query_end++;
        if (ref_winnow_end < ref_hash.minimizers.size() && ref_hash.minimizers[ref_winnow_end].second < ref_end + 1) { // extend reference
            winnow.add_to_reference(ref_hash.minimizers[ref_winnow_end].first);
            ref_winnow_end++;
        } 
        ref_end++;

        if (winnow.valid_jaccard()) {
            prev_jaccard = double(winnow.jaccard) / winnow.query.size();
            prev_jaccard_p = winnow.jaccard;
            best_ref_end = ref_end;
            best_query_end = query_end;
            recover = 0;
        } else {
            recover++;
        }
    }

E:  edist = filter(query_hash.seq, query_start, 
        best_query_end - query_start + 1, ref_hash.seq, 
        ref_start, best_ref_end - ref_start + 1
    );
    if (edist.first.first < 0) {
        return;
    }

    hits.push_back({
        query_start, best_query_end, ref_start, best_ref_end, 
        j2md(prev_jaccard), j2md(0), 0, // fix these 
        make_pair(0, prev_jaccard_p), edist.first
    });
    
    auto a = INTERVAL(query_start, best_query_end);
    auto b = INTERVAL(ref_start, best_ref_end);
    TREE += make_pair(a, SUBTREE_t({b, {make_pair(a, b)}}));
}

vector<Hit> search (int query_start, 
                    const Hash &ref_hash, 
                    const Hash &query_hash, 
                    bool allow_overlaps) 
{
    if (!qgram_p.size()) {
        qgram_p = vector<int>(1 << (2 * 8), 0);
        qgram_r = vector<int>(1 << (2 * 8), 0);
        aho = new AHOAutomata();
    }
 
    int st = query_hash.find_minimizers(query_start);
    if (st == query_hash.minimizers.size()) 
        return vector<Hit>();
    int mi = st;
    
    auto pf = TREE.find(query_start);

    // Iterate through all unique hashes in query[query_start: query_start + MIN_READ_SIZE]
    // `candidates` is a list of positions in the reference which match the query hashes
    SlidingMap<hash_t> init_winnow;
    vector<int> candidates;
    for (; mi < query_hash.minimizers.size() && query_hash.minimizers[mi].second - query_start <= MIN_READ_SIZE; mi++) { 
        auto &h = query_hash.minimizers[mi].first;
        if (init_winnow.find(h) != init_winnow.end())
            continue;
        init_winnow.add_to_query(h);

        if (h.first) // If it is N hash, ignore it
            continue;
        
        // Is this hash in the reference? Is this hash high-call hash as well?
        auto ptr = ref_hash.index.find(h);
        if (ptr == ref_hash.index.end() || ptr->second.size() >= ref_hash.threshold) 
            continue;

        // Iterate through positions in the reference with this hash
        for (auto &pos: ptr->second) {
            // Make sure to have at least 1 kb spacing if reference = query
            if (!allow_overlaps && pos < query_start + MIN_READ_SIZE)
                continue;
            // if (!check_overlap(pf, query_start, pos)) {
            //     continue;
            // }
            candidates.push_back(pos);
        }
    }
    sort(candidates.begin(), candidates.end());
    
    if (!init_winnow.size()) { // TODO min sketch size limit for Ns
        return vector<Hit>();
    }
    
    int M = ceil(init_winnow.size() * tau());
    // for (auto e: candidates) eprn("{:n}", e);
    // eprn(">> search: {} candidates to evaluate (M={}, init_winnow_query={}, tau={})", candidates.size(), M, init_winnow_query.size(), tau());

    // Find all locations in the `candidates` so that
    // MIN_READ_SIZE read covers at least M (= s * tau) hashes
    auto T = vector<pair<int, int>>{{-1, -1}};
    for (int i = 0; i <= (int)candidates.size() - M; i++) {
        int j = i + (M - 1);
        if (candidates[j] - candidates[i] < MIN_READ_SIZE) {
            int x = max(0, candidates[j] - MIN_READ_SIZE + 1), 
                y = candidates[i] + 1;
            if (x >= T.back().second)
                T.push_back({x, y});
            T.back().second = y;
        }
    }

    vector<Hit> hits;
    for (int t = 1; t < T.size(); t++) {
        if (!allow_overlaps)
            T[t].first = max(T[t].first, query_start + MIN_READ_SIZE);
        if (T[t].first > T[t].second) 
            continue;

        eprn(">> extend query:{:n} vs ref:{:n}...{:n}", query_start, T[t].first, T[t].second);
    
        TOTAL_ATTEMPTED++;
        auto winnow = init_winnow;

        // auto pf = TREE.find(query_start);
        // if (!check_overlap(pf, query_start, min_ref)) {
        //     INTERVAL_FAILED ++;
        //     return;
        // }
        
        int ref_start = T[t].first, 
            ref_end   = T[t].first + MIN_READ_SIZE;
        int ref_winnow_start = ref_hash.find_minimizers(ref_start);
        assert(ref_winnow_start < ref_hash.minimizers.size());
        
        // winnow is W(query) ; extend it to W(query) | W(ref) and mark elements in W(query) & W(ref)
        int ref_winnow_end = ref_winnow_start;
        for (; ref_winnow_end < ref_hash.minimizers.size() && ref_hash.minimizers[ref_winnow_end].second < ref_end; ref_winnow_end++) {
            winnow.add_to_reference(ref_hash.minimizers[ref_winnow_end].first);
        }
        winnow.rewind();

        // eprn(">> extend needed_jaccard:{} jaccard:{} init_start:{} init_jaccard={}", 
        //     tau * winnow_query.size(),
        //     jaccard, 
        //     seed_ref_start, seed_jaccard);
        // Roll until we find best inital match
        bool found = false;
        while (ref_start < T[t].second) {
            if (ref_winnow_start < ref_hash.minimizers.size() && ref_hash.minimizers[ref_winnow_start].second < ref_start + 1) {
                winnow.remove_from_reference(ref_hash.minimizers[ref_winnow_start].first);
                ref_winnow_start++;    
            }
            if (ref_winnow_end < ref_hash.minimizers.size() && ref_hash.minimizers[ref_winnow_end].second < ref_end + 1) {
                winnow.add_to_reference(ref_hash.minimizers[ref_winnow_end].first);
                ref_winnow_end++;
            }
            if (winnow.valid_jaccard()) {
                extend(winnow, 
                    query_hash, query_start, query_start + MIN_READ_SIZE, st, mi,
                    ref_hash, ref_start, ref_end, ref_winnow_start, ref_winnow_end, 
                    hits, allow_overlaps
                );
                found = true;
            }
            ref_start++, ref_end++;
        }

        if (!found) JACCARD_FAILED ++;
    }
    TREE -= INTERVAL(0, query_start);

    return parse_hits(hits);
}

