/// 786

#include <bits/stdc++.h>
#include "common.h"
#include "search.h"
using namespace std;

#include <boost/icl/interval_map.hpp>
#include <boost/icl/interval_set.hpp>

/// http://www.biorxiv.org/content/biorxiv/early/2017/03/24/103812.full.pdf

const int MAX_MATCH = 1 * 1000 * 1000; // 1 MB at most!

typedef boost::icl::discrete_interval<int> INTERVAL;
boost::icl::interval_map<int, boost::icl::interval_set<int>> TREE;


inline bool in_map_n(const auto &m, const hash_t &k)
{ 
    return !k.first && m.find(k) != m.end();
}

inline pair<int, bool> in_map_nc(const auto &m, const hash_t &k, int ref_pos)
{ 
    if (k.first) return {ref_pos, false};
    auto p = m.find(k);
    if (p == m.end()) return {ref_pos, false};
    if (ref_pos == -1) return {p->second.first, true};
    return {ref_pos, true};
}


void refine(int p, int idx_p, int idx_q, // query start and W(query) range
            int x, int y, // window to check for the initial MIN_READ_SIZE match
            int initial_match_length,
            //map<hash_t, char> 
            auto L0, // W(query)
            const Hash &ref_hash, 
            const Hash &query_hash,
            vector<Hit> &hits, // Output
            bool allow_overlaps=false)
{
    // chr22   16050000        16051868        chr22   48680215        48682083
    double tau = ::tau();

    // prn("*** {} to {}/{} ***", p, x, y);
    auto pf = TREE.find(p);
    boost::icl::interval_set<int>::iterator pfp;
    if (pf != TREE.end() && (pfp = pf->second.find(x)) != pf->second.end()) { // check overlap!
        return;
    }

    int i = x, 
        j = x + initial_match_length;

    // L0 is |W(query)| and L0[h] is {location of h in the query, false}
    // L  is |W(reference) + W(query)| 
    // L[h] is {location of h in the reference, true} iff h is in |W(reference) & W(query)|
    SlidingMap<hash_t, pair<int, bool>> L; // matches: ref---query!
    for (auto &p: L0) L.store[p.first] = {-1, false};
    int idx_i = ref_hash.find_minimizers(i);
    if (idx_i == ref_hash.minimizers.size()) {
        return;
    }
    int idx_j = idx_i;
    while (idx_j < ref_hash.minimizers.size() && ref_hash.minimizers[idx_j].second < j) {
        // do not set 1 to N hashes
        L.store[ref_hash.minimizers[idx_j].first] = in_map_nc(L0, ref_hash.minimizers[idx_j].first, ref_hash.minimizers[idx_j].second);
        idx_j++;
    }

    L.boundary = next(L.store.begin(), L0.size() - 1);
    int jaccard = L.boundary->second.second;
    for (auto it = L.store.begin(); it != L.boundary; it++)
        jaccard += it->second.second;
    // prn("ja={}|||{}", jaccard, L0.size());
    // Extend it up to MIN_READ_SIZE (initial match)
    int seed_i = -1, seed_jaccard;
    if (jaccard >= tau * L0.size())
        seed_i = i, seed_jaccard = jaccard;
    if (seed_i != -1) while (i < y) {
        if (idx_i < ref_hash.minimizers.size() && ref_hash.minimizers[idx_i].second < i + 1) {
            jaccard += L.remove(ref_hash.minimizers[idx_i].first);
            idx_i++;        
        }
        if (idx_j < ref_hash.minimizers.size() && ref_hash.minimizers[idx_j].second < j + 1) {
            jaccard += L.add(ref_hash.minimizers[idx_j].first, in_map_nc(L0, ref_hash.minimizers[idx_j].first, ref_hash.minimizers[idx_j].second));
            idx_j++;
        }
        if (jaccard >= tau * L0.size() && jaccard >= seed_jaccard) {
            seed_i = i, seed_jaccard = jaccard;
            break;
        }
        i++, j++;
    }
    if (seed_i == -1)
        return;
    // prn("ja={} se_i={} se_ja={}", jaccard, seed_i, seed_jaccard);

    // Keep extending it as much as we can

    // L has all minimizers now
    // we are mapping query[p, q] to ref[i, j]
    int q = p + initial_match_length;
    // idx_i/idx_j points to the first/last minimizer of REF_i,j
    // idx_p/idx_q points to the first/last minimizer of QUERY_p,q
    double prev_jaccard = double(jaccard) / L0.size();
    double init_jaccard = prev_jaccard;
    int last_n1 = 0, 
        last_n2 = 0;

    const int RECOVER_BP = 250; // We allow 250bp extra extend just in case!
    int recover = 0;

    int best_q = q;
    int best_j = j;
    char reason = 0;
    vector<pair<int, int>> best_matches;
    for (auto &p: L.store) if (p.second.second) {
        assert(L0.find(p.first) != L0.end());
        best_matches.push_back({p.second.first, L0[p.first].first}); // ref-query
    }
    // prn("init best match: {}", best_matches.size());


    while (q < query_hash.seq.size() && j < ref_hash.seq.size()) {
        // prn(">> p,q={},{} ({}) i,j={},{} ({}) jaccard={} s={} score={}",
        //     p,q,q-p,i,j,j-i,jaccard,L0.size(), j2md(double(jaccard)/L0.size()));
    
    // - disallow overlaps or too long matches
        int max_match = MAX_MATCH;
        if (!allow_overlaps) 
            max_match = min(max_match, int((1.0 / ERROR_RATE + .5) * abs(p - i)));
        if (max(q - p, j - i) > max_match) {
            // hits.push_back({p, q - 1, i, j - 1,  j2md(prev_jaccard), 0});
            reason = 0;
            break;
        }

    // - try extending to the right
        // extend query
        if (query_hash.minimizers[idx_q].second < q + 1) {
            auto h = query_hash.minimizers[idx_q].first;
            if (!in_map(L0, h)) {
                L0[h] = {query_hash.minimizers[idx_q].second, false};
                // check is h in L now
                jaccard += L.add(h, in_map_nc(L.store, h, -1));
                L.boundary++; // |W(A)| increases
                jaccard += L.boundary->second.second;
            }
            idx_q++;
        }
        q++;
        // extend reference
        if (idx_j < ref_hash.minimizers.size() && ref_hash.minimizers[idx_j].second < j + 1) { 
            auto h = ref_hash.minimizers[idx_j];
            jaccard += L.add(h.first, in_map_nc(L0, h.first, h.second));
            idx_j++;
        } 
        j++;
        if (jaccard >= tau * L0.size()) {
            prev_jaccard = double(jaccard) / L0.size();
            best_j = j;
            best_q = q;
            recover = 0;
            best_matches.clear();
            for (auto &p: L.store) if (p.second.second) {
                assert(L0.find(p.first) != L0.end());
                best_matches.push_back({p.second.first, L0[p.first].first}); // ref-query
            }
            continue;
        } 

        // 1. keep extending some times...
        if (recover < RECOVER_BP) {
            recover++;
            continue;
        }

        reason = 2;
        break;
    }

    // prn("bmatch # {}", best_matches.size());
    hits.push_back({p, best_q, i, best_j, j2md(prev_jaccard), j2md(init_jaccard), reason, best_matches});
    TREE += make_pair(INTERVAL(p, best_q), 
            boost::icl::interval_set<int>({INTERVAL(i, best_j)}));
    if (allow_overlaps) 
        TREE += make_pair(INTERVAL(i, best_j),
            boost::icl::interval_set<int>({INTERVAL(p, best_q)}));
}

vector<Hit> search (int query_start, 
                    const Hash &ref_hash, 
                    const Hash &query_hash, 
                    bool allow_overlaps) 
{
    map<hash_t, pair<int, bool>> L0;
    int st = query_hash.find_minimizers(query_start);
    if (st == query_hash.minimizers.size())
        return vector<Hit>();
    int mi = st;
    //prn("{}\n", mi);
    
    auto pf = TREE.find(query_start);
    boost::icl::interval_set<int>::iterator pfp;

    // Iterate through all unique hashes in query[query_start: query_start + MIN_READ_SIZE]
    // `candidates` is a list of positions in the reference which match the query hashes
    vector<int> candidates;
    for (; mi < query_hash.minimizers.size() && query_hash.minimizers[mi].second - query_start <= MIN_READ_SIZE; mi++) { 
        auto &m = query_hash.minimizers[mi];
        if (mi != st && m.first == query_hash.minimizers[mi - 1].first) 
            continue;
        
        L0[m.first] = {m.second, false};
        // If it is N hash, ignore it
        if (m.first.first)
            continue;
        // Is this hash in the reference? Is this hash high-call hash as well?
        auto ptr = ref_hash.index.find(m.first);
        if (ptr == ref_hash.index.end() || ptr->second.size() >= ref_hash.threshold) 
            continue;
        // Iterate through positions in the reference with this hash
        for (auto &pos: ptr->second) {
            // Make sure to have at least 1 kb spacing if reference = query
            // TODO: handle reverse complements
            if (!allow_overlaps && pos < query_start + 2 * MIN_READ_SIZE)
                continue;

            int initial_match_length = MIN_READ_SIZE;
            if (pf != TREE.end() && (pfp = pf->second.find(pos)) != pf->second.end()) // check overlap!
            {
                continue; // should be extended on the right later on!

                // prn("{} - {}", pf->first.lower(), pf->first.upper());
                // prn("{} < {}", pfp->lower(), pfp->upper());

                // int sA = pf->first.lower(), eA = pf->first.upper();
                // int sB = pfp->lower(),      eB = pfp->upper();

                // they need to have AT LEAST
                // int endA = eA - query_start; 
                // assert(endA > 0);
                // int endB = eB - pos; 
                // assert(endB > 0);
                // initial_match_length = min(endA, endB) + MIN_READ_SIZE;
                
                // if (initial_match_length > MAX_MATCH) continue;

                // if (abs(abs(sA - query_start) - abs(sB - pos)) < MIN_READ_SIZE) {
                //     continue;
                // } else {
                //     initial_match_length = min(abs(eA - query_start), abs(eB - pos));
                // }
            }

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
    vector<pair<int, int>> T = {{-1, -1}};
    for (int i = 0; i <= (int)candidates.size() - M; i++) {
        int j = i + (M - 1);

        // int newM = ceil(candidates[i].second * tau());
        // int j = i + (newM - 1);

        if (candidates[j] - candidates[i] < MIN_READ_SIZE) {
            int x = max(0, candidates[j] - MIN_READ_SIZE + 1), 
                y = candidates[i] + 1;
            if (x >= T.back().second)
                T.push_back({x, y});
            T.back().second = y;
        }
    }

    vector<Hit> hits, hits_real;
    // for (int t = 1; t < T.size(); t++) 
    //     prn("{} {}-{}", query_start,  T[t].first, T[t].second);
    // Extend each candidate!
    for (int t = 1; t < T.size(); t++) {
        refine(
            query_start, st, mi, 
            T[t].first, T[t].second, 
            MIN_READ_SIZE,
            L0, ref_hash, query_hash, hits, allow_overlaps
        );
    }
    sort(hits.begin(), hits.end(), [](const auto &a, const auto &b) {
        if (a.q != b.q) return a.q > b.q; // .p are all equal
        return tie(a.i, a.j) < tie(a.i, a.j);   // larger intervals at the top! then by other ints
    });
    for (auto &h: hits) { // brute force, but maybe easier than doing full interval thing
            auto &pp=h;
            // prn("{}\t{}\t{}\t{}\t{}\t{}\t\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tBM;{}",
            //     ">>>>>", pp.p, pp.q, 
            //     "<<<<<", pp.i, pp.j,
            //     pp.id, 
            //     "+", "+",
            //     // Optional fields
            //     int(pp.break_criteria),
            //     pp.q - pp.p, pp.j - pp.i,
            //     pp.init_id,
            //     pp.matches.size()
            // );

        bool add = true;
        for (int j = hits_real.size() - 1; j >= 0; j--) {
            auto &ph = hits_real[j];
            if (h.i >= ph.i && h.j <= ph.j) { // if full match
                add = false;
                break;
            }
        }
        if (add) hits_real.push_back(h);
    }
    // hits = hits_real;
    // for (auto &h: hits_real) {
    //     TREE += make_pair(INTERVAL(h.p, h.q), 
    //         boost::icl::interval_set<int>({INTERVAL(h.i, h.j)}));
    // }
    reverse(hits_real.begin(), hits_real.end()); // smaller to larger
    TREE -= INTERVAL(0, query_start);

    return hits_real;
}

