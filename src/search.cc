/// 786

#include <bits/stdc++.h>
#include "common.h"
#include "search.h"
using namespace std;

/// http://www.biorxiv.org/content/biorxiv/early/2017/03/24/103812.full.pdf

void refine(int p, int idx_p, int idx_q, // query start and W(query) range
            int x, int y, // window to check for the initial MIN_READ_SIZE match
            map<uint32_t, char> L0, // W(query)
            const Hash &ref_hash, 
            const Hash &query_hash,
            vector<Hit> &hits, // Output
            bool allow_overlaps=true);

vector<Hit> search (int query_start, 
                    const Hash &ref_hash, 
                    const Hash &query_hash, 
                    bool allow_overlaps) 
{
    vector<int> candidates;

    map<uint32_t, char> L0;
    uint64_t prev_hash = 1LL << 63;
    int st = query_hash.find_minimizers(query_start);
    int mi = st;
    for (; query_hash.minimizers[mi].second - query_start <= MIN_READ_SIZE; mi++) { 
        /// ^^ TODO count initial kmer overlap

        auto &m = query_hash.minimizers[mi];
        if (m.first == prev_hash) continue;
        prev_hash = m.first;
        
        L0[m.first] = 0;
        auto ptr = ref_hash.index.find(m.first);
        if (ptr == ref_hash.index.end() || ptr->second.size() >= ref_hash.threshold) 
            continue;
        for (auto &pos: ptr->second) {
            // Make sure to have at least 1 kb spacing if reference = query
            // TODO: handle reverse complements
            if (!allow_overlaps && (!(pos < query_start - MIN_READ_SIZE || pos > query_start + 2 * MIN_READ_SIZE)))
                continue;
            candidates.push_back(pos);
        }
    }
    sort(candidates.begin(), candidates.end());
    
    int M = ceil(L0.size() * tau());
    // prn("{} candidates to evaluate (M={})", candidates.size(), M);

    vector<pair<int, int>> T = {{-1, -1}};
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
    for (auto &t: T) {
        refine(query_start, st, mi, t.first, t.second, L0, ref_hash, query_hash, hits, allow_overlaps);
    }
    return hits;
}

void refine(int p, int idx_p, int idx_q, // query start and W(query) range
            int x, int y, // window to check for the initial MIN_READ_SIZE match
            map<uint32_t, char> L0, // W(query)
            const Hash &ref_hash, 
            const Hash &query_hash,
            vector<Hit> &hits, // Output
            bool allow_overlaps)
{
    double tau = ::tau();

    int i = x, 
        j = x + MIN_READ_SIZE;

    SlidingMap<uint32_t, char> L(L0);
    int idx_i = ref_hash.find_minimizers(i);
    int idx_j = idx_i;
    while (ref_hash.minimizers[idx_j].second < j) {
        L.store[ref_hash.minimizers[idx_j].first] = in_map(L0, ref_hash.minimizers[idx_j].first);
        idx_j++;
    }

    L.boundary = next(L.store.begin(), L0.size() - 1);
    int jaccard = L.boundary->second;
    for (auto it = L.store.begin(); it != L.boundary; it++)
        jaccard += it->second;
    
// Extend it up to MIN_READ_SIZE (initial match)
    int seed_i = -1, seed_jaccard;
    if (jaccard >= tau * L0.size())
        seed_i = i, seed_jaccard = jaccard;
    if (seed_i != -1) while (i < y) {
        if (ref_hash.minimizers[idx_i].second < i + 1) {
            jaccard += L.remove(ref_hash.minimizers[idx_i].first);
            idx_i++;        
        }
        if (ref_hash.minimizers[idx_j].second < j + 1) {
            jaccard += L.add(ref_hash.minimizers[idx_j].first, in_map(L0, ref_hash.minimizers[idx_j].first));
            idx_j++;
        }
        if (jaccard >= tau * L0.size() && jaccard >= seed_jaccard) {
            seed_i = i, seed_jaccard = jaccard;
            break;
        }
        i++, j++;
    }

// Keep extending it as much as we can

    // L has all minimizers now
    // we are mapping query[p, q] to ref[i, j]
    int q = p + MIN_READ_SIZE;
    // idx_i/idx_j points to the first/last minimizer of REF_i,j
    // idx_p/idx_q points to the first/last minimizer of QUERY_p,q
    int prev_jaccard = jaccard;
    int last_n1 = 0, 
        last_n2 = 0;
    while (true) {
        // prn(">> p,q={},{} ({}) i,j={},{} ({}) jaccard={} s={} score={}",p,q,q-p,i,j,j-i,jaccard,L0.size(), j2md(double(jaccard)/L0.size(), KMER_SIZE));
    
    // - if references are the same, disallow overlaps
    //   TODO: handle reverse complements!
        if (!allow_overlaps && ((i < q && p < j) || (q - p > 1000000))) {
            hits.push_back({p, q - 1, i, j - 1,  j2md(double(prev_jaccard) / L0.size()), 0});
            break;
        }

    // - if we have more than 10 Ns in either query or, terminate
        if (query_hash.seq[q - 1] == 'N') last_n1++; else last_n1 = 0;
        if (ref_hash.seq[j - 1] == 'N') last_n2++; else last_n2 = 0;
        if (last_n1 > 10 || last_n2 > 10) {
            hits.push_back({p, q - 1, i, j - 1, j2md(double(prev_jaccard) / L0.size()), 1});
            break;
        }

    // - try extending to the right
        // extend query
        if (query_hash.minimizers[idx_q].second < q + 1) {
            auto h = query_hash.minimizers[idx_q].first;
            if (!in_map(L0, h)) {
                L0[h] = 0;
                jaccard += L.add(h, in_map(L.store, h));
                L.boundary++; // |W(A)| increases
                jaccard += L.boundary->second;
            }
            idx_q++;
        }
        q++;
        // extend reference
        if (ref_hash.minimizers[idx_j].second < j + 1) { 
            auto h = ref_hash.minimizers[idx_j].first;
            jaccard += L.add(h, in_map(L0, h));
            idx_j++;
        } 
        j++;
        if (jaccard >= tau * L0.size()) {
            prev_jaccard = jaccard;
            continue;
        }

    // - if jaccard is not enough, remember previous case and try shrinking it from the left
    //   TODO: better criteria for this...
        hits.push_back({p, q - 1, i, j - 1, j2md(double(prev_jaccard) / L0.size()), 2});
        // shrink query
        if (query_hash.minimizers[idx_p].second < p + 1) {
            auto h = query_hash.minimizers[idx_p].first;
            L0.erase(h);
            jaccard += L.remove(h);
            jaccard -= L.boundary->second;
            L.boundary--;
            idx_p++;
        }
        p++;
        // shrink reference
        if (ref_hash.minimizers[idx_i].second < i + 1) {
            auto h = ref_hash.minimizers[idx_i].first;
            jaccard += L.remove(h);
            idx_i++;
        }
        i++;
        if (jaccard >= tau * L0.size()) {
            prev_jaccard = jaccard;
            continue;
        }

    // - if all fails, stop searching
        break;
    }
}

