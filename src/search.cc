/// 786

#include <bits/stdc++.h>
#include "common.h"
#include "search.h"
using namespace std;

#include <boost/icl/interval_map.hpp>
#include <boost/icl/interval_set.hpp>
#include <edlib.h>

/// http://www.biorxiv.org/content/biorxiv/early/2017/03/24/103812.full.pdf

/*
504a505
> chr1  40020441    40021507    chr1    40021508    40022596
512a514
> chr1  55002548    55003748    chr1    241083029   241084244
*/

const int MAX_MATCH = 1 * 500 * 1000; // 1 MB at most!

typedef boost::icl::discrete_interval<int> INTERVAL;
boost::icl::interval_map<int, boost::icl::interval_set<int>> TREE;

int QGRAM_NORMAL_FAILED = 0;
int QGRAM_SPACED_FAILED = 0;
int EDIT_FAILED = 0;

void printAlignment(const char* query, const char* target,
                    const unsigned char* alignment, const int alignmentLength,
                    const int position, const EdlibAlignMode modeCode) 
{
    int tIdx = -1;
    int qIdx = -1;
    if (modeCode == EDLIB_MODE_HW) {
        tIdx = position;
        for (int i = 0; i < alignmentLength; i++) {
            if (alignment[i] != EDLIB_EDOP_INSERT)
                tIdx--;
        }
    }
    for (int start = 0; start < alignmentLength; start += 50) {
        // target
        eprnn("T: ");
        int startTIdx;
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            if (alignment[j] == EDLIB_EDOP_INSERT)
                eprnn("-");
            else
                eprnn("{}", target[++tIdx]);
            if (j == start)
                startTIdx = tIdx;
        }
        eprnn(" ({} - {})\n", max(startTIdx, 0), tIdx);

        // match / mismatch
        eprnn("   ");
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            eprnn(alignment[j] == EDLIB_EDOP_MATCH ? "|" : " ");
        }
        eprnn("\n");

        // query
        eprnn("Q: ");
        int startQIdx = qIdx;
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            if (alignment[j] == EDLIB_EDOP_DELETE)
                eprnn("-");
            else
                eprnn("{}", query[++qIdx]);
            if (j == start)
                startQIdx = qIdx;
        }
        eprnn(" ({} - {})\n\n", max(startQIdx, 0), qIdx);
    }
}

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

inline int min_qgram(int l, int q) {
    return l * (1 - MAX_GAP_ERROR - q * MAX_EDIT_ERROR) - (GAP_FREQUENCY * l + 1) * (q - 1);
}

inline int min_qgram_lo(int l, int q) {
    return l * (1 - MAX_GAP_ERROR - q * (MAX_EDIT_ERROR/2)) - (GAP_FREQUENCY * l + 1) * (q - 1);
}


vector<int> qgram_p, qgram_r;
auto filter(const string &q, int q_pos, int q_len, const string &r, int r_pos, int r_len) 
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
    //     // int QG = get<1>(patterns[pi]);
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
    //         // eprn("wohooo {}", pattern);
    //         return make_pair(-1, -1);
    //     }
    // }

    // vector<int> qgram_p(QSZ, 0);
    // vector<int> qgram_i(QSZ, 0);

    // const int MASK = QSZ - 1;
    // for (int qi = q_pos, qgram = 0; qi < q_pos + q_len; qi++) {
    //     qgram = ((qgram << 2) | qdna(q[qi])) & MASK;
    //     if (qi - q_pos >= QG) qgram_p[qgram] += 1;
    // }
    // for (int qi = r_pos, qgram = 0; qi < r_pos + r_len; qi++) {
    //     qgram = ((qgram << 2) | qdna(r[qi])) & MASK;
    //     if (qi - r_pos >= QG) qgram_i[qgram] += 1;
    // }
    // int dist = 0;
    // for (int qi = 0; qi < QSZ; qi++) 
    //     dist += min(qgram_p[qi], qgram_i[qi]);

    // int maxlen = max(q_len, r_len);
    // int maxedit = maxlen / 4 + 10;
    // if (dist < maxlen - QSZ * maxedit - (QSZ - 1))
    //     return make_pair(-1, -1);


    int edist = 0;
    // auto result = edlibAlign(
    //     q.c_str() + q_pos, q_len,
    //     r.c_str() + r_pos, r_len,
    //     edlibNewAlignConfig(maxlen * .40, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, 0, 0)
    // );
    // auto edist = result.editDistance;
    // edlibFreeAlignResult(result);
    // if (edist == -1) {
    //     return -1;
    // }

    return make_pair(0, 0);
    // if (result.editDistance == -1 > DL/4 && dist >= DL - QG + 1 - QG * (DL / 4)) {
    //     eprn("\n--- {}: ed {} qg {}/{} jac {}", DL, result.editDistance, 
    //         dist, DL - QG + 1 - QG * (DL / 4), prev_jaccard_p);
    //     // char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
    //     printAlignment(q.c_str() + q_pos, ref_hash.seq.c_str() + r_pos, result.alignment, 
    //         result.alignmentLength, *result.endLocations, EDLIB_MODE_NW);
    //     exit(0);
    // }
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
// > chr22 21033606 21034790    chr22   21038016    21039210
// chr22   21032636    21038019    chr22   21417024    21411669
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
    SlidingMap<hash_t, bool> L(L0); // matches: ref---query!
    // for (auto &p: L0) L.store[p.first] = {-1, false};
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

    L.boundary = next(L.store.begin(), L0.size() - 1);
    int jaccard = L.boundary->second;
    for (auto it = L.store.begin(); it != L.boundary; it++)
        jaccard += it->second;
    // prn("ja={}|||{}", jaccard, L0.size());
    // Extend it up to MIN_READ_SIZE (initial match)
    int seed_i = -1, seed_jaccard = 0;
    if (jaccard >= tau * L0.size())
        seed_i = i, seed_jaccard = jaccard;
    // prn("[{} {}] ja={} se_i={} se_ja={}", x, y, jaccard, seed_i, seed_jaccard);
    //if (seed_i != -1) 
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
    if (seed_i == -1)
        return;

    // Keep extending it as much as we can

    // L has all minimizers now
    // we are mapping query[p, q] to ref[i, j]
    int q = p + initial_match_length;
    // idx_i/idx_j points to the first/last minimizer of REF_i,j
    // idx_p/idx_q points to the first/last minimizer of QUERY_p,q
    double prev_jaccard = double(jaccard) / L0.size();
    double init_jaccard = prev_jaccard;
    int prev_jaccard_p = jaccard;
    int init_jaccard_p = jaccard;
    int last_n1 = 0, 
        last_n2 = 0;

    const int RECOVER_BP = 250; // We allow 250bp extra extend just in case!
    int recover = 0;

    int best_q = q;
    int best_j = j;
    char reason = 0;
    // vector<pair<int, int>> best_matches;
    // for (auto &p: L.store) if (p.second.second) {
    //     assert(L0.find(p.first) != L0.end());
    //     best_matches.push_back({p.second.first, L0[p.first].first}); // ref-query
    // }
    // prn("init best match: {}", best_matches.size());

    while (q < query_hash.seq.size() && j < ref_hash.seq.size()) {
        // prn(">> p,q={},{} ({}) i,j={},{} ({}) jaccard={} s={} score={}",
        //     p,q,q-p,i,j,j-i,jaccard,L0.size(), j2md(double(jaccard)/L0.size()));
    
        // - disallow overlaps or too long matches
        int max_match = MAX_MATCH;
        if (!allow_overlaps) 
            max_match = min(max_match, int((1.0 / MAX_GAP_ERROR + .5) * abs(p - i)));
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
                L0[h] = 0; // {query_hash.minimizers[idx_q].second, false};
                // check is h in L now
                jaccard += L.add(h, in_map_n(L.store, h));
                L.boundary++; // |W(A)| increases
                jaccard += L.boundary->second;
            }
            idx_q++;
        }
        q++;
        // extend reference
        if (idx_j < ref_hash.minimizers.size() && ref_hash.minimizers[idx_j].second < j + 1) { 
            auto &h = ref_hash.minimizers[idx_j];
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
            // best_matches.clear();
            // for (auto &p: L.store) if (p.second.second) {
            //     assert(L0.find(p.first) != L0.end());
            //     best_matches.push_back({p.second.first, L0[p.first].first}); // ref-query
            // }
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

    auto edist = filter(query_hash.seq, p, best_q - p + 1, ref_hash.seq, i, best_j - i + 1);
    if (edist.first < 0)
        return;

    // prn("bmatch # {}", best_matches.size());
    hits.push_back({p, best_q, i, best_j, j2md(prev_jaccard), j2md(init_jaccard), reason, make_pair(init_jaccard_p, prev_jaccard_p), edist.second});
    TREE += make_pair(INTERVAL(p, best_q), 
            boost::icl::interval_set<int>({INTERVAL(i, best_j)}));
    // if (allow_overlaps) 
    //     TREE += make_pair(INTERVAL(i, best_j),
    //         boost::icl::interval_set<int>({INTERVAL(p, best_q)}));
}

vector<Hit> search (int query_start, 
                    const Hash &ref_hash, 
                    const Hash &query_hash, 
                    bool allow_overlaps) 
{

    if (!qgram_p.size()) {
        qgram_p = vector<int>(1<<(2*8),0);
        qgram_r = vector<int>(1<<(2*8),0);
    }

    map<hash_t, bool> L0;
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
        
        L0[m.first] = false;
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

    // eprn("{} candidates to refine", T.size());

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
            //     pp.id, cat 
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

