/// 786
/// Fast winnowing algorithm: 
//  https://people.cs.uct.ac.za/~ksmith/articles/sliding_window_minimum.html

/******************************************************************************/

#include <vector>
#include <string>
#include <deque>
#include <algorithm>

#include "common.h"
#include "hash.h"

using namespace std;

/******************************************************************************/

ostream& operator<<(ostream& os, const Hash& dt)
{  
	os << int(dt.status) << '.' << std::hex << dt.hash;
	return os;  
}  

bool operator<(const Hash &x, const Hash &y) 
{
	return tie(x.status, x.hash) < tie(y.status, y.hash);
}

bool operator==(const Hash &x, const Hash &y) 
{
	return tie(x.status, x.hash) == tie(y.status, y.hash);
}

bool operator!=(const Hash &x, const Hash &y) 
{
	return tie(x.status, x.hash) != tie(y.status, y.hash);
}

bool operator<=(const Hash &x, const Hash &y) 
{
	return (x < y) || (x == y);
}

bool operator==(const Minimizer &x, const Minimizer &y) 
{
	return tie(x.loc, x.hash) == tie(y.loc, y.hash);
}

bool operator<(const Minimizer &x, const Minimizer &y) 
{
	return tie(x.loc, x.hash) < tie(y.loc, y.hash);
}

/******************************************************************************/


vector<Minimizer> get_minimizers(const string &s, int kmer_size, 
	const int window_size, bool separate_lowercase)
{
	// int SKIP = 3;
	// const int kmer_span = (kmer_size / SKIP) * (SKIP + 1);
	// assert(kmer_span < s.size());
	// assert(kmer_span < 32);
	// const uint64_t MASK  = (1ull << (2 * kmer_span)) - 1;
	// const uint32_t MASKN = (1ul << kmer_span) - 1;
	// const uint64_t UMASK = (1ull << (2 * SKIP)) - 1;

	vector<Minimizer> minimizers;
	minimizers.reserve((2 * s.size()) / window_size);
	deque<Minimizer> window;

	// uint32_t last_n = 0;
	// uint32_t last_u = 0;
	// for (int i = 0; i < kmer_span; i++) {
	// 	h = ((h << 2) | hash_dna(s[i])) & MASK;
	// 	last_n = ((last_n << 1) | (s[i] == 'N')) & MASKN;
	// 	last_u = ((last_u << 1) | (s[i] != 'N' && isupper(s[i]))) & MASKN;
	// }

	// // window contains k-mer *starting positions* (i.e. the last k-mer's end might go outside of the window)
	// for (int i = 0; i < s.size() - kmer_span; i++) {
	// 	uint32_t hash1 = h >> (2 * (kmer_span - kmer_size));
	// 	uint32_t hash2 = 0;
	// 	hash2 = (hash2 << (2 * SKIP)) | ((h >> (2 * (3 * (SKIP + 1) + 1))) & UMASK);
	// 	hash2 = (hash2 << (2 * SKIP)) | ((h >> (2 * (2 * (SKIP + 1) + 1))) & UMASK);
	// 	hash2 = (hash2 << (2 * SKIP)) | ((h >> (2 * (1 * (SKIP + 1) + 1))) & UMASK);
	// 	hash2 = (hash2 << (2 * SKIP)) | ((h >> (2 * (0 * (SKIP + 1) + 1))) & UMASK);
	// 	// bool has_u = (last_u >= i);
	// 	bool has_n_1 = (last_n >> (kmer_span - kmer_size));
	// 	bool has_n_2 = (last_n & 0xEEEE);
	// 	bool has_u_1 = (last_u >> (kmer_span - kmer_size));
	// 	bool has_u_2 = (last_u & 0xEEEE);

	// 	// uint32_t y = 0;
	// 	// for (int j = 0; j < kmer_size; j++)
	// 		// y = ((y << 2) | hash_dna(s[i + j + j/3])) & ((1ull << (2 * kmer_size)) - 1);
	// 	// assert(y==hash2);
	// 	// hash2 = y;

	// 	// for (int j = 0; j < kmer_span; j++) eprnn("{}", s[i+j]);  eprnn(" ");
	// 	// for (int j = 0; j < kmer_span; j++) eprnn("{}", (last_u>>(kmer_span-j-1))&1);  eprnn(" ");
	// 	// for (int j = 0; j < kmer_size; j++) eprnn("{}", s[i+j]);  eprnn(" ");
	// 	// for (int j = 0; j < kmer_size; j++) eprnn("{}", s[i+j+j/3]); eprn(" --> {} {}", has_u_1, has_u_2);
		
	// 	uint32_t hc = hash1;
	// 	bool has_n = has_n_1; 
	// 	bool has_u = has_u_1;
	// 	if (tie(has_n_2, has_u_1, hash2) < tie(has_n_1, has_u_2, hash1)) {
	// 		// if (has_u_1 != has_u_2)
	// 		// 	eprn("{:x} {} {} -> {:x} {} {}", hc, has_n, has_u, hash2, has_n_2, has_u_2);
	// 		hc = hash2, has_n = has_n_2, has_u = has_u_2;
	// 	}

	// 	h = ((h << 2) | hash_dna(s[i + kmer_span])) & MASK;
	// 	last_n = ((last_n << 1) | (s[i + kmer_span] == 'N')) & MASKN;
	// 	last_u = ((last_u << 1) | (s[i + kmer_span] != 'N' && isupper(s[i + kmer_span]))) & MASKN;

	// 	Hash hh { hc, has_n 
	// 		? Hash::Status::HAS_N 
	// 		: (has_u ? Hash::Status::HAS_UPPERCASE : Hash::Status::ALL_LOWERCASE)
	// 	};   
	// 	// eprn("{}", hh);
	// 	if (!separate_lowercase && hh.status == Hash::Status::ALL_LOWERCASE) {
	// 		hh.status = Hash::Status::HAS_UPPERCASE;
	// 	}
	// 	while (!window.empty() && !(window.back().hash < hh)) {
	// 		window.pop_back();
	// 	}
	// 	while (!window.empty() && window.back().loc < i - window_size) {
	// 		window.pop_front();
	// 	}
	// 	window.push_back({hh, i});

	// 	if (i < window_size) 
	// 		continue;
	// 	if (!minimizers.size() || !(window.front() == minimizers.back())) {
	// 		minimizers.push_back(window.front());
	// 	}
	// }
	// // exit(0);
	// return minimizers;

	const uint32_t MASK  = (1ull << (2 * kmer_size)) - 1;
	uint32_t h = 0;
	int last_n = - kmer_size - window_size;
	int last_u = last_n;
	for (int i = 0; i < s.size()-3; i++) {
		if (s[i] == 'N') {
			last_n = i;
		} else if (isupper(s[i])) {
			last_u = i;
		}

		h = ((h << 2) | hash_dna(s[i])) & MASK; 

		if (i < kmer_size) 
			continue;

		// uint32_t cur_h = ((h << 2) | hash_dna(s[i+2]));
		// cur_h = ((cur_h << 2) | hash_dna(s[i+3]));
		Hash hh { h, last_n >= (i - kmer_size + 1) 
			? Hash::Status::HAS_N 
			: (last_u >= (i - kmer_size + 1)) // || isupper(s[i+2]) || isupper(s[i+3])
				? Hash::Status::HAS_UPPERCASE 
				: Hash::Status::ALL_LOWERCASE
		};   
		if (!separate_lowercase && hh.status == Hash::Status::ALL_LOWERCASE) {
			hh.status = Hash::Status::HAS_UPPERCASE;
		}
		while (!window.empty() && !(window.back().hash < hh)) {
			window.pop_back();
		}
		while (!window.empty() && window.back().loc < (i - kmer_size + 1) - window_size) {
			window.pop_front();
		}
		window.push_back({hh, i - kmer_size + 1});

		if (i - kmer_size + 1 < window_size) 
			continue;
		if (!minimizers.size() || !(window.front() == minimizers.back())) {
			minimizers.push_back(window.front());
		}
	}

	// sort(minimizers.begin(), minimizers.end());
	return minimizers;
}

/******************************************************************************/

Sequence::Sequence(const string &name, const string &seq, bool is_rc):
	name(name), seq(seq), is_rc(is_rc)
{
	if (is_rc) {
		this->seq = rc(seq);
	}
}

/******************************************************************************/

Index::Index(shared_ptr<Sequence> seq, int kmer_size, int window_size, bool separate_lowercase): 
	seq(seq), kmer_size(kmer_size), window_size(window_size) 
{
	// eprn("Hashing {} bps", s.size());

	assert(kmer_size <= 16);
	minimizers = get_minimizers(seq->seq, kmer_size, window_size, separate_lowercase);
	
	for (auto &i: minimizers) {
		index[i.hash].push_back(i.loc);
	}

	int ignore = (minimizers.size() * 0.001) / 100.0;
	
	map<int, int> hist;
	for (auto &i: index) {
		hist[i.second.size()] += 1;
	}
	int sum = 0;
	threshold = 1 << 31;
	int j = 0;
	for (auto i = hist.rbegin(); i != hist.rend(); i++, j++) {
		sum += i->second;
		if (sum <= ignore) { 
			threshold = i->first;
		} else { 
			break;
		}
	}
	// eprn("Ignoring top {} hits", j);
}

int Index::find_minimizers(int p) const
{
	int lo = 0, hi = minimizers.size() - 1, mid;
	while (lo <= hi) {
		mid = lo + (hi - lo) / 2;
		if (minimizers[mid].loc >= p && (!mid || minimizers[mid - 1].loc < p))
			break;
		if (minimizers[mid].loc < p) {
			lo = mid + 1;
		} else {
			hi = mid;
		}
	}
	assert(minimizers[mid].loc >= p || mid == minimizers.size() - 1); 
	assert(!mid || minimizers[mid-1].loc < p);
	if (minimizers[mid].loc < p) {
		mid++; // the last one--- no solution
	}
	return mid;
}
