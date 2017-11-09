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

const int MAX_MATCH = 500 * 1000; // 50kb at most
const int RIGHT_ALLOWANCE = 750; /// TODO more mathy formulation

int64_t TOTAL_ATTEMPTED = 0;
int64_t JACCARD_FAILED = 0;
int64_t QGRAM_NORMAL_FAILED = 0;
int64_t OTHER_FAILED = 0;
int64_t CORE_FAILED = 0;
int64_t INTERVAL_FAILED = 0;

inline int min_qgram(int l, int q) {
    return l * (1 - MAX_GAP_ERROR - q * MAX_EDIT_ERROR) - (GAP_FREQUENCY * l + 1) * (q - 1);
}

inline bool hits_cmp(const Hit &a, const Hit &b) {
    if (a.q != b.q) return a.q > b.q; // .query_start are all equal
    return make_pair(a.i, a.j) < make_pair(a.i, a.j);   // larger intervals at the top! then by other ints
}

inline bool check_overlap (TREE_t::const_iterator &pf, int pf_pos, int pfp_pos) {
    SUBTREE_t::const_iterator pfp;
    if (pf != TREE.end() && (pfp = pf->second.find(pfp_pos)) != pf->second.end()) { // check overlap!
        for (set<pair<INTERVAL, INTERVAL> >::iterator it = pfp->second.begin(); it != pfp->second.end(); it++) {
            int sA = it->first.lower(),  eA = it->first.upper();
            int sB = it->second.lower(), eB = it->second.upper();
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

inline void add_to_tree(INTERVAL a, INTERVAL b) {
    set<pair<INTERVAL, INTERVAL> > s;
    s.insert(make_pair(a, b));
    boost::icl::interval_map<int, set<pair<INTERVAL, INTERVAL> > > iss;
    iss += make_pair(b, s);
    TREE += make_pair(a, iss);
}

AHOAutomata *aho = NULL;
vector<int> qgram_p, qgram_r;

pair<pair<int, int>, pair<int, int>> filter(const string &q, int q_pos, int q_len, const string &r, int r_pos, int r_len) 
{
    int q_up = 0; for (int i = 0; i < q_len; i++) q_up += isupper(q[q_pos + i]); 
    int r_up = 0; for (int i = 0; i < r_len; i++) r_up += isupper(r[r_pos + i]);

    if (q_up < 500 || r_up < 500) {
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
    for (map<int, int>::iterator it = hits.begin(); it != hits.end(); it++)
        if (it->second == 3) common++;

    double boundary = (1/2) * (q_len / 50.0);
    if (common < boundary) {
        CORE_FAILED++;
        return make_pair(make_pair(-1, -1), make_pair(-1, -1));
    }

    return make_pair(make_pair(dist, common), make_pair(q_up, r_up));
}

inline bool match_winnow_query(multimap<hash_t, bool> &m, const hash_t &k, multimap<hash_t, bool>::iterator &it) { 
    if (k.first) return false;

    it = m.lower_bound(k);
    for (; it != m.end() && it->first == k; it++) {
        if (!it->second) return true;
    }
    return false;
}

void refine(int query_start, int query_winnow_end, // query start and W(query) range
            int min_ref, int max_ref, // window to check for the initial MIN_READ_SIZE match
            //map<hash_t, char> 
            multimap<hash_t, bool> winnow_query, // W(query)
            const Hash &ref_hash, 
            const Hash &query_hash,
            vector<Hit> &hits, // Output
            bool allow_overlaps=false)
{
    double tau = ::tau();

    // eprn(">> eval {} vs {}..{}", query_start, min_ref, max_ref);

    TREE_t::const_iterator pf = TREE.find(query_start);
    if (!check_overlap(pf, query_start, min_ref)) {
        TOTAL_ATTEMPTED ++;
        INTERVAL_FAILED ++;
        return;
    }

    multimap<hash_t, bool>::iterator tmpit;
    
    int ref_start = min_ref, 
        ref_end = min_ref + MIN_READ_SIZE;
    // winnow_query is |W(query)| and winnow_query[h] is {location of h in the query, false}
    // winnow_union  is |W(reference) + W(query)| 
    // winnow_union[h] is {location of h in the reference, position} iff h is in |W(reference) & W(query)|
    SlidingMap<hash_t, bool> winnow_union(winnow_query); // matches: ref -> query!
    int ref_winnow_start = ref_hash.find_minimizers(ref_start);
    if (ref_winnow_start == ref_hash.minimizers.size()) {
        TOTAL_ATTEMPTED ++;
        JACCARD_FAILED ++;
        return;
    }
    int ref_winnow_end = ref_winnow_start;
    // winnow_union is W(query) as winnow_query; 
    // now extend it to W(query)|W(ref) and mark elements in W(query)&W(ref)

    // prnn("winnow query init");
    // for (auto &it : winnow_query) {
    //     prnn("{}{}:{} ", "*N"[it.first.first],it.first.second, it.second);
    // } prnn("\n");

    while (ref_winnow_end < ref_hash.minimizers.size() && ref_hash.minimizers[ref_winnow_end].second < ref_end) {
        // do not set 1 to N hashes
        // prn("checking {}{}", "*N"[ref_hash.minimizers[ref_winnow_end].first.first],ref_hash.minimizers[ref_winnow_end].first.second);
        bool is_in = match_winnow_query(winnow_query, ref_hash.minimizers[ref_winnow_end].first, tmpit);
        auto it = winnow_union.store.insert(make_pair(ref_hash.minimizers[ref_winnow_end].first, is_in));
        if (is_in) tmpit->second = true; // it is marked as NOT TOUCHED [fix: pointer maybe?]
        ref_winnow_end++;
    }
    // eprnn("init sha={}\n", shared_init);

    winnow_union.boundary = winnow_union.store.begin();
    std::advance(winnow_union.boundary, winnow_query.size() - 1);
    
    int jaccard = winnow_union.boundary->second;
    for (multimap<hash_t, bool>::iterator it = winnow_union.store.begin(); it != winnow_union.boundary; it++)
        jaccard += it->second;

    // Extend it up to MIN_READ_SIZE (initial match)
    int seed_ref_start = -1, seed_jaccard = 0;
    if (jaccard >= tau * winnow_query.size())
        seed_ref_start = ref_start, seed_jaccard = jaccard;

    int ref_uppercase = 0;
    int query_uppercase = 0;

    // prn("[{} {}] ja={} se_i={} se_ja={}", min_ref, max_ref, jaccard, seed_ref_start, seed_jaccard);
    while (ref_start < max_ref) {
        if (ref_winnow_start < ref_hash.minimizers.size() && ref_hash.minimizers[ref_winnow_start].second < ref_start + 1) {
            // removing reference hash 
            auto &h = ref_hash.minimizers[ref_winnow_start].first;
            // LOWER BOUND: remove matching hash "just in case"
            auto itu = winnow_union.store.lower_bound(h); 
            bool in_union = false;
            if (itu != winnow_union.store.end() && itu->first == h) {
                auto itu2 = itu;
                for (; itu2 != winnow_union.store.end() && itu2->first == h; itu2++)
                    if (itu2->second) {
                        itu = itu2; 
                        break;
                    }

                // itu is "the hash"
                // delink it from winnow_query
                tmpit = winnow_query.lower_bound(h);
                for (; tmpit != winnow_query.end() && tmpit->first == h; tmpit++)
                    if (tmpit->second) {
                        tmpit->second = false;
                        break;
                    }

                if (h < winnow_union.boundary->first) {
                    winnow_union.boundary++;
                    jaccard += winnow_union.boundary->second - itu->second;
                } else if (h == winnow_union.boundary->first) {
                    itu2 = itu;
                    for (; itu2 != winnow_union.store.end() && itu2->first == h; itu2++)
                        if (itu2 == winnow_union.boundary) {
                            winnow_union.boundary++;
                            jaccard += winnow_union.boundary->second - itu->second;
                            break;
                        }
                }
                winnow_union.store.erase(itu);
            }            
            ref_winnow_start++;    
        }
        if (ref_winnow_end < ref_hash.minimizers.size() && ref_hash.minimizers[ref_winnow_end].second < ref_end + 1) {
            bool is_in = match_winnow_query(winnow_query, ref_hash.minimizers[ref_winnow_end].first, tmpit);
            jaccard += winnow_union.add(ref_hash.minimizers[ref_winnow_end].first, is_in);
            if (is_in) tmpit->second = true;
            ref_winnow_end++;
        }
        if (jaccard >= tau * winnow_query.size() && jaccard >= seed_jaccard) {
            seed_ref_start = ref_start, seed_jaccard = jaccard;
            break;
        }
        ref_start++, ref_end++;
    }
    if (seed_ref_start == -1) {
        TOTAL_ATTEMPTED ++;
        JACCARD_FAILED ++;
        return;
    }

    // Now keep extending it as much as we can
    // winnow_union has all minimizers now
    // we are mapping query[query_start, q] to ref[ref_start, ref_end]
    int query_end = query_start + MIN_READ_SIZE;
    // ref_winnow_start/ref_winnow_end points to the first/last minimizer of ref
    // query_winnow_start/query_winnow_end points to the first/last minimizer of query
    double prev_jaccard = double(jaccard) / winnow_query.size();
    double init_jaccard = prev_jaccard;
    int prev_jaccard_p = jaccard;
    int init_jaccard_p = jaccard;
    
    const int RECOVER_BP = 250; // We allow 250bp extra extend just in case!
    int recover = 0;

    int best_query_end = query_end;
    int best_ref_end = ref_end;
    char reason = 0;
    
    auto edist = filter(query_hash.seq, query_start, best_query_end - query_start + 1, 
        ref_hash.seq, ref_start, best_ref_end - ref_start + 1);
    if (edist.first.first < 0) {
        TOTAL_ATTEMPTED ++;
        return;
    }

    int qlow = best_query_end - query_start - edist.second.first;
    int rlow = best_ref_end - ref_start - edist.second.second;

    while (query_end < query_hash.seq.size() && ref_end < ref_hash.seq.size()) {    
        if (qlow > (query_end-query_start)/4 || rlow > (ref_end-ref_start)/4) {
            break;
        }

        // - disallow overlaps or too long matches
        int max_match = MAX_MATCH;
        if (!allow_overlaps) 
            max_match = min(max_match, int((1.0 / MAX_GAP_ERROR + .5) * abs(query_start - ref_start)));
        if (max(query_end - query_start, ref_end - ref_start) > max_match) {
            reason = 0;
            break;
        }

        // N detection
        int n_count = 0;
        while ((query_hash.seq[query_end] == 'N') || (ref_hash.seq[ref_end] == 'N')) {
            n_count++, query_end++, ref_end++;
            if (n_count > 100) goto end;
        }

        if (max(ref_uppercase, query_uppercase) < (ref_end - ref_start) / 10) {
            goto end;
        }


        // - try extending to the right
        if (query_hash.minimizers[query_winnow_end].second < query_end + 1) { // extend query
            hash_t h = query_hash.minimizers[query_winnow_end].first;
            // check is h in winnow_union now
            auto it = winnow_query.insert(make_pair(h, 0)); // {query_hash.minimizers[query_winnow_end].second, false};
            bool is_in = match_winnow_query(winnow_union.store, h, tmpit);
            if (is_in) it->second = tmpit->second = true; // mark hash in both winnows as matched!
            jaccard += winnow_union.add(h, is_in);
            winnow_union.boundary++; // |W(A)| increases
            jaccard += winnow_union.boundary->second;
            query_winnow_end++;
        }
        if (islower(query_hash.seq[query_end])) qlow++;
        query_end++;
        if (ref_winnow_end < ref_hash.minimizers.size() && ref_hash.minimizers[ref_winnow_end].second < ref_end + 1) { // extend reference
            const minimizer_t &h = ref_hash.minimizers[ref_winnow_end];
            bool is_in = match_winnow_query(winnow_query, h.first, tmpit);
            jaccard += winnow_union.add(h.first, is_in);
            if (is_in) tmpit->second = true;
            ref_winnow_end++;
        } 
        if (islower(ref_hash.seq[ref_end])) rlow++;
        ref_end++;
        if (jaccard >= tau * winnow_query.size()) {
            prev_jaccard = double(jaccard) / winnow_query.size();
            prev_jaccard_p = jaccard;
            best_ref_end = ref_end;
            best_query_end = query_end;
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

    end:
    edist = filter(query_hash.seq, query_start, best_query_end - query_start + 1, ref_hash.seq, ref_start, best_ref_end - ref_start + 1);
    if (edist.first.first < 0) {
        TOTAL_ATTEMPTED ++;
        return;
    }

    TOTAL_ATTEMPTED ++;

    hits.push_back(Hit(query_start, best_query_end, ref_start, best_ref_end, j2md(prev_jaccard), j2md(init_jaccard), reason, make_pair(init_jaccard_p, prev_jaccard_p), edist.first));
    add_to_tree(INTERVAL(query_start, best_query_end), INTERVAL(ref_start, best_ref_end));
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
 
    multimap<hash_t, bool> winnow_query;
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
        
        winnow_query.insert(make_pair(m.first, 0));
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
            if (!check_overlap(pf, query_start, pos)) {
                TOTAL_ATTEMPTED ++;
                INTERVAL_FAILED ++;
                continue;
            }
            candidates.push_back(pos);
        }
    }
    sort(candidates.begin(), candidates.end());
    
    if (winnow_query.size() < 1) { // TODO min sketch size limit for Ns
        TOTAL_ATTEMPTED ++;
        JACCARD_FAILED ++;
        return vector<Hit>();
    }
    int M = ceil(winnow_query.size() * tau());
    // eprn("{} candidates to evaluate (M={}, winnow_query={}, tau={})", candidates.size(), M, winnow_query.size(), tau());

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
            query_start, mi, // st, mi, 
            T[t].first, T[t].second, 
            winnow_query, ref_hash, query_hash, hits, allow_overlaps
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

