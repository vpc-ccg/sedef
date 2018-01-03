/// 786
/// Fast winnowing algorithm: 
//  https://people.cs.uct.ac.za/~ksmith/articles/sliding_window_minimum.html

/******************************************************************************/

#include <vector>
#include <string>
#include <deque>

#include "common.h"
#include "hash.h"

using namespace std;

/******************************************************************************/

ostream& operator<<(ostream& os, const hash_t& dt)
{  
	os << int(dt.first) << '.' << std::hex << (unsigned long long)dt.second;  
	return os;  
}  

/******************************************************************************/

vector<minimizer_t> get_minimizers(const string &s, const int kmer_size, const int window_size)
{
	const uint32_t MASK = (1 << (2 * kmer_size)) - 1;

	vector<minimizer_t> minimizers;
	minimizers.reserve((2 * s.size()) / window_size);
	deque<minimizer_t> window;
	int last_n = - kmer_size - window_size;
	int last_u = last_n;
	// window here is defined as list of WINDOW_SIZE 
	// k-mer *starting positions* (i.e. last k-mer goes outside of the window)
	for (uint32_t i = 0, h = 0; i < s.size(); i++) {
		if (s[i] == 'N') {
			last_n = i;
		} else if (isupper(s[i])) {
			last_u = i;
		}

		h = ((h << 2) | hash_dna(s[i])) & MASK; 
		if (i < kmer_size) 
			continue;

		char status = 1 - char(last_u >= (i - kmer_size + 1)); // 1 for no-uppercase, 0 for uppercase
		if (last_n >= (i - kmer_size + 1)) 
			status = 2;
		hash_t hh = make_pair(status, h);   
		while (!window.empty() && (window.back().first >= hh)) {
			window.pop_back();
		}
		while (!window.empty() && window.back().second < (i - kmer_size + 1) - window_size) {
			window.pop_front();
		}
		window.push_back(make_pair(hh, i - kmer_size + 1));

		if (i - kmer_size + 1 < window_size) 
			continue;
		if (!minimizers.size() || window.front() != minimizers.back()) {
			minimizers.push_back(window.front());
		}
	}
	return minimizers;
}

/******************************************************************************/

Hash::Hash(const string &s, int kmer_size, int window_size): 
	seq(s), kmer_size(kmer_size), window_size(window_size) 
{
	// eprn("Hashing {} bps", s.size());

	assert(kmer_size <= 16);
	minimizers = get_minimizers(s, kmer_size, window_size);
	
	for (auto &i: minimizers) {
		index[i.first].push_back(i.second);
	}

	int ignore = (minimizers.size() * 0.001) / 100.0;
	
	map<int, int> hist;
	for (auto &i: index) {
		hist[i.second.size()] += 1;
	}
	int sum = 0;
	threshold = 1 << 31;
	for (auto i = hist.rbegin(); i != hist.rend(); i++) {
		sum += i->second;
		if (sum <= ignore) { 
			threshold = i->first;
		} else { 
			break;
		}
	}
	// eprn("Index cut-off threshold: {}", threshold);
}

int Hash::find_minimizers(int p) const
{
	int lo = 0, hi = minimizers.size() - 1, mid;
	while (lo <= hi) {
		mid = lo + (hi - lo) / 2;
		if (minimizers[mid].second >= p && (!mid || minimizers[mid - 1].second < p))
			break;
		if (minimizers[mid].second < p) {
			lo = mid + 1;
		} else {
			hi = mid;
		}
	}
	assert(minimizers[mid].second >= p || mid == minimizers.size() - 1); 
	assert(!mid || minimizers[mid-1].second < p);
	if (minimizers[mid].second < p) {
		mid++; // the last one--- no solution
	}
	return mid;
}
