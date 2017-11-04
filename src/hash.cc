/// 786

#include <vector>
#include <string>
#include <deque>
#include "common.h"
#include "hash.h"
using namespace std;

vector<minimizer_t> get_minimizers(const string &s)
{
    vector<minimizer_t> minimizers;
    minimizers.reserve((2 * s.size()) / WINDOW_SIZE);
    deque<minimizer_t> window;
    uint32_t h = 0;
    const uint32_t MASK = (1 << (2 * KMER_SIZE)) - 1;
    int last_n = - KMER_SIZE - WINDOW_SIZE;
    for (int i = 0; i < s.size(); i++) {
        if (s[i] == 'N') last_n = i;

        h = (h << 2) | qdna(s[i]); 
        h &= MASK;
        if (i < KMER_SIZE) continue;

        hash_t hh = make_pair(bool(last_n >= (i - KMER_SIZE + 1) - WINDOW_SIZE), h);   
        while (!window.empty() && (window.back().first >= hh))
            window.pop_back();
        while (!window.empty() && window.back().second <= (i - KMER_SIZE + 1) - WINDOW_SIZE)
            window.pop_front();
        // Do not add minimizers with N
        window.push_back(make_pair(hh, i - KMER_SIZE + 1));

        if (i - KMER_SIZE + 1 < WINDOW_SIZE) continue;
        if (!minimizers.size() || window.front() != minimizers.back()) {
            minimizers.push_back(window.front());
        }
    }
    return minimizers;
}

Hash::Hash (const string &s): seq(s)
{
    eprn("Hashing {} bps", s.size());
    minimizers = get_minimizers(s);
    
    for (int i = 0; i < minimizers.size(); i++) {
        index[minimizers[i].first].push_back(minimizers[i].second);
    }

    int ignore = (minimizers.size() * 0.001) / 100.0;
    
    map<int, int> hist;
    for (map<hash_t, list<int>, MapCompare>::iterator i = index.begin(); i != index.end(); i++)
        hist[i->second.size()] += 1;
    int sum = 0;
    threshold = 1 << 31;
    for (map<int, int>::reverse_iterator i = hist.rbegin(); i != hist.rend(); i++) {
        sum += i->second;
        if (sum <= ignore) threshold = i->first;
        else break;
    }
    eprn("Index cut-off threshold: {}", threshold);
}

int Hash::find_minimizers(int p) const
{
    int lo = 0, hi = minimizers.size() - 1, mid;
    while (lo <= hi) {
        mid = lo + (hi - lo) / 2;
        if (minimizers[mid].second >= p && (!mid || minimizers[mid - 1].second < p))
            break;
        if (minimizers[mid].second < p) lo = mid + 1;
        else hi = mid;
    }
    assert(minimizers[mid].second >= p || mid == minimizers.size() - 1); // if last one, well, what can we do...
    assert(!mid || minimizers[mid-1].second < p);
    if (minimizers[mid].second < p) mid++; // zadnji, nema rjesenja
    return mid;
}