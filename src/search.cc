/// 786
/// Based on http://www.biorxiv.org/content/biorxiv/early/2017/03/24/103812.full.pdf

#include <list>
#include <queue>
#include <vector>
#include <string>
#include <cmath>
#include <queue>
#include "common.h"
#include "filter.h"
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
int64_t INTERVAL_FAILED = 0;

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

auto parse_hits(vector<Hit> &hits)
{
    // COUNT SQUEEZED!
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

void extend(SlidingMap &winnow, SlidingMap2 &winnow2,
    const Hash &query_hash, int query_start, int query_end,
    int query_winnow_start, int query_winnow_end,
    const Hash &ref_hash, int ref_start, int ref_end,
    int ref_winnow_start, int ref_winnow_end,
    vector<Hit> &hits,
    bool allow_overlaps) 
{
    eprn(">> extend query:{:n} vs ref:{:n}", query_start, ref_start);
            
    if (!filter(query_hash.seq, query_start, query_end - query_start + 1, 
            ref_hash.seq, ref_start, ref_end - ref_start + 1)) 
    {
        return;
    }

    int qlow = 0; //query_end - query_start - edist.second.first;
    int rlow = 0; //ref_end - ref_start - edist.second.second;

    int best_query_end = query_end;
    int best_ref_end = ref_end;
    int j = winnow.jaccard();
    double prev_jaccard = double(j) / winnow.query_size;
    int prev_jaccard_p = j;
    for (int recover = 0; query_end < query_hash.seq.size() && ref_end < ref_hash.seq.size() && recover < RECOVER_BP; ) {   
        //break;
        // if (qlow > (query_end-query_start)/4 || rlow > (ref_end-ref_start)/4) { // TODO FIX
        //     break;
        // }
        // if (islower(query_hash.seq[query_end])) qlow++;
        // if (islower(ref_hash.seq[ref_end])) rlow++;

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
        while (query_end < query_hash.seq.size() && ref_end < ref_hash.seq.size() 
                && ((query_hash.seq[query_end] == 'N') || (ref_hash.seq[ref_end] == 'N'))) 
        {
            n_count++, query_end++, ref_end++;
            if (n_count > 100) goto E;
        }

        // - try extending to the right
        if (query_winnow_end < query_hash.minimizers.size() 
                && query_hash.minimizers[query_winnow_end].second < query_end + 1) 
        { // extend query
            winnow.add_to_query(query_hash.minimizers[query_winnow_end].first);
            winnow2.add_to_query(query_hash.minimizers[query_winnow_end].first);
            assert(winnow.jaccard() == winnow2.jaccard());
            query_winnow_end++;
        }
        if (ref_winnow_end < ref_hash.minimizers.size() 
                && ref_hash.minimizers[ref_winnow_end].second < ref_end + 1) 
        { // extend reference
            winnow.add_to_reference(ref_hash.minimizers[ref_winnow_end].first);
            winnow2.add_to_reference(ref_hash.minimizers[ref_winnow_end].first);
            assert(winnow.jaccard() == winnow2.jaccard());
            ref_winnow_end++;
        } 
        if ((j = winnow.jaccard()) >= 0) {
            prev_jaccard = double(j) / winnow.query_size;
            prev_jaccard_p = j;
            best_ref_end = ref_end;
            best_query_end = query_end;
            recover = 0;
        } else {
            recover++;
        }

        query_end++, ref_end++;
    }

E:  if (!filter(query_hash.seq, query_start, best_query_end - query_start + 1, 
            ref_hash.seq, ref_start, best_ref_end - ref_start + 1))
    {
        return;
    }
    
    hits.push_back({
        query_start, best_query_end, 
        ref_start, best_ref_end, 
        j2md(prev_jaccard), prev_jaccard_p
    });
    
    auto a = INTERVAL(query_start, best_query_end);
    auto b = INTERVAL(ref_start, best_ref_end);
    // TREE += make_pair(a, SUBTREE_t({b, {make_pair(a, b)}}));
}

vector<Hit> search (int query_start, 
                    const Hash &ref_hash, 
                    const Hash &query_hash, 
                    bool allow_overlaps) 
{ 
    int st = query_hash.find_minimizers(query_start);
    if (st == query_hash.minimizers.size()) 
        return vector<Hit>();
    int mi = st;
    
    auto pf = TREE.find(query_start);

    // Iterate through all unique hashes in query[query_start: query_start + MIN_READ_SIZE]
    // `candidates` is a list of positions in the reference which match the query hashes
    SlidingMap init_winnow;
    SlidingMap2 init_winnow2;
    set<int> candidates_prel;
    for (; mi < query_hash.minimizers.size() && query_hash.minimizers[mi].second - query_start <= MIN_READ_SIZE; mi++) { 
        auto &h = query_hash.minimizers[mi].first;
        
        init_winnow.add_to_query(h);
        init_winnow2.add_to_query(h);
        
        if (h.first) // If it is N hash, ignore it
            continue;
        
        auto ptr = ref_hash.index.find(h);
        if (ptr == ref_hash.index.end() || ptr->second.size() >= ref_hash.threshold) 
            continue;
        for (auto &pos: ptr->second) {
            // Make sure to have at least 1 kb spacing if reference = query
            if (!allow_overlaps && pos < query_start + MIN_READ_SIZE)
                continue;
            // if (!check_overlap(pf, query_start, pos)) {
            //     continue;
            // }
            candidates_prel.insert(pos);
        }
    }
    if (!init_winnow.query_size) { // TODO min sketch size limit for Ns
        return vector<Hit>();
    }
    int M = ceil(init_winnow.query_size * tau());
    eprn("==M {}", M);
    // int M = init_winnow.size() * tau();
    
    // for (auto e: candidates) eprn("{:n}", e);
    // eprn(">> search: {} candidates to evaluate (M={}, init_winnow_query={}, tau={})", candidates.size(), M, init_winnow_query.size(), tau());

    // Find all locations in the `candidates` so that
    // MIN_READ_SIZE read covers at least M (= s * tau) hashes
    auto T = vector<pair<int, int>>{{-1, -1}};
    //sort(candidates.begin(), candidates.end());
    vector<int> candidates(candidates_prel.begin(), candidates_prel.end());
    for (int i = 0; i <= (int)candidates.size() - M; i++) {
        int j = i + (M - 1);
        if (candidates[j] - candidates[i] < MIN_READ_SIZE) {
            int x = max(0, candidates[j] - MIN_READ_SIZE + 1), 
                y = candidates[i] + 1;
            //if (x >= T.back().second)
            T.push_back({x, y});
            //T.back().second = y; // 
        }
    }

    vector<Hit> hits;
    for (int t = 1; t < T.size(); t++) {
        if (!allow_overlaps)
            T[t].first = max(T[t].first, query_start + MIN_READ_SIZE);
        if (T[t].first > T[t].second) 
            continue;

        // eprn(">> extend query:{:n} vs ref:{:n}...{:n}", query_start, T[t].first, T[t].second);
    
        TOTAL_ATTEMPTED++;
        auto winnow = init_winnow;
        auto winnow2 = init_winnow2;

        // auto pf = TREE.find(query_start);
        // if (!check_overlap(pf, query_start, min_ref)) {
        //     INTERVAL_FAILED ++;
        //     return;
        // }
        
        int ref_start = T[t].first, 
            ref_end   = T[t].first + MIN_READ_SIZE;
        int ref_winnow_start = ref_hash.find_minimizers(ref_start);
        assert(ref_winnow_start < ref_hash.minimizers.size());

        auto px = [&](int x) {
            string s1;
            if (x==1){
                for (auto e: winnow) s1 += fmt::format("{}:{} ", (int)e.first.first.first, e.first.first.second); 
            } else {
                map<hash_t, int> wwq, wwr, wwp;
                for (auto x: winnow2.query) wwq[x]++, wwp[x]=0;
                for (auto x: winnow2.ref) wwr[x]++, wwp[x]=0;
                for (auto &x: wwp) x.second = max(wwq[x.first], wwr[x.first]);
                for (auto &e: wwp) for(int i=0;i<e.second;i++) s1 += fmt::format("{}:{} ", 
                    (int)e.first.first, e.first.second); 
            }
            return s1;
        };

        winnow.rewind();
        winnow2.rewind();
        assert(winnow.jaccard() == winnow2.jaccard());
        assert(px(1)==px(2));
        
        // winnow is W(query) ; extend it to W(query) | W(ref) and mark elements in W(query) & W(ref)
        int ref_winnow_end = ref_winnow_start;
        for (; ref_winnow_end < ref_hash.minimizers.size() && ref_hash.minimizers[ref_winnow_end].second < ref_end; ref_winnow_end++) {
            winnow.add_to_reference(ref_hash.minimizers[ref_winnow_end].first);
            winnow2.add_to_reference(ref_hash.minimizers[ref_winnow_end].first);
            assert(winnow.jaccard() == winnow2.jaccard());

            eprn("ai {}:{}", (int)ref_hash.minimizers[ref_winnow_end].first.first, ref_hash.minimizers[ref_winnow_end].first.second);
            if (px(1)!=px(2)) eprn("{}\n{}", px(1), px(2));
            assert(px(1)==px(2));
        }

        // eprn(">> extend needed_jaccard:{} jaccard:{} init_start:{} init_jaccard={}", 
        //     tau * winnow_query.size(),
        //     jaccard, 
        //     seed_ref_start, seed_jaccard);
        // Roll until we find best inital match
        int best_jaccard = winnow.jaccard();
        int best_ref_start = ref_start, best_ref_end = ref_end;
        int best_ref_winnow_start = ref_winnow_start, best_ref_winnow_end = ref_winnow_end;
        assert(px(1)==px(2));
        while (ref_start < T[t].second) {
            if (ref_winnow_start < ref_hash.minimizers.size() && ref_hash.minimizers[ref_winnow_start].second < ref_start + 1) {
                winnow.remove_from_reference(ref_hash.minimizers[ref_winnow_start].first);
                winnow2.remove_from_reference(ref_hash.minimizers[ref_winnow_start].first);

                eprn("\nrm {}:{}", (int)ref_hash.minimizers[ref_winnow_start].first.first, ref_hash.minimizers[ref_winnow_start].first.second);
                eprn("-- {} vs true {}", winnow.jaccard() , winnow2.jaccard());
                if (px(1)!=px(2)) eprn("{}\n{}", px(1), px(2));
                assert(px(1)==px(2));
                assert(winnow.jaccard() == winnow2.jaccard());
                ref_winnow_start++;    
            }
            if (ref_winnow_end < ref_hash.minimizers.size() && ref_hash.minimizers[ref_winnow_end].second < ref_end + 1) {
                winnow.add_to_reference(ref_hash.minimizers[ref_winnow_end].first);
                winnow2.add_to_reference(ref_hash.minimizers[ref_winnow_end].first);
                eprn("\nad {}:{}", (int)ref_hash.minimizers[ref_winnow_end].first.first, ref_hash.minimizers[ref_winnow_end].first.second);
                eprn("-- {} vs true {}", winnow.jaccard() , winnow2.jaccard());
                if (px(1)!=px(2)) eprn("{}\n{}", px(1), px(2));
                assert(px(1)==px(2));
                assert(winnow.jaccard() == winnow2.jaccard());
                ref_winnow_end++;
            }
            int j;
            if ((j = winnow.jaccard()) > best_jaccard) {
                best_jaccard = j;
                best_ref_start = ref_start;
                best_ref_end = ref_end;
                best_ref_winnow_start = ref_winnow_start;
                best_ref_winnow_end = ref_winnow_end;
            }
            ref_start++, ref_end++;
        }

        assert(winnow.jaccard() == winnow2.jaccard());
        if (best_jaccard < 0) JACCARD_FAILED++;
        else extend(winnow, winnow2,
            query_hash, query_start, query_start + MIN_READ_SIZE, st, mi,
            ref_hash, best_ref_start, best_ref_end, best_ref_winnow_start, best_ref_winnow_end, 
            hits, allow_overlaps
        );
    }
   
    TREE -= INTERVAL(0, query_start);
    return parse_hits(hits);
}

