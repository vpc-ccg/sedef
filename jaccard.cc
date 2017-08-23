#include <bits/stdc++.h>
// #define FMT_HEADER_ONLY
#include <fmt/format.h>
using namespace std;
using namespace fmt;

constexpr auto dna_lookup_init() {
	array<char, 128> values = {};
	get<'C'>(values) = 1; get<'c'>(values) = 1;
	get<'G'>(values) = 2; get<'g'>(values) = 2;
	get<'T'>(values) = 3; get<'t'>(values) = 3;
	return values;
}
constexpr auto dna_lookup = dna_lookup_init();
inline char qdna(char c) {
	return dna_lookup[c];
}

typedef pair<uint32_t, int> minimizer_t;
typedef unordered_map<uint32_t, list<int>> hash_t;

const int KMER_SIZE = 16;
static_assert(KMER_SIZE <= 16, "k-mer space is 32-bit");

const int WINDOW_SIZE = 16; // <-- Needs to be changed
const int MIN_ALIGNMENT = 1000;
const int SKETCH_SIZE = 2 * MIN_ALIGNMENT / WINDOW_SIZE;
const int MIN_READ_SIZE = 1000;

chrono::time_point<chrono::high_resolution_clock> t_start, t_end;

auto get_minimizers(string s) 
{
	vector<minimizer_t> minimizers;
	minimizers.reserve((2 * s.size()) / WINDOW_SIZE);
	deque<minimizer_t> window;
	uint32_t h = 0;
	for (int i = 0; i < s.size(); i++) {
		h = (h << 2) | qdna(s[i]); 
		if (i < KMER_SIZE) continue;
		
		while (!window.empty() && window.back().first >= h)
			window.pop_back();
		while (!window.empty() && window.back().second <= (i - KMER_SIZE + 1) - WINDOW_SIZE)
			window.pop_front();
		window.push_back({h, i - KMER_SIZE + 1});

		if (i - KMER_SIZE + 1 < WINDOW_SIZE) continue;
		if (!minimizers.size() || window.front() != minimizers.back()) {
			minimizers.push_back(window.front());
		}
	}
	return minimizers;
}

// k-mer -> (pos)
auto get_hash (string s) 
{
	auto minimizers = get_minimizers(s);
	hash_t hash;
	for (auto &x: minimizers) {
		hash[x.first].push_back(x.second);
	}
	return make_pair(minimizers, hash);
}

inline double tau(double error = .75)
{
	return 1 / (2 * exp(KMER_SIZE * error) - 1);
}

auto get_minimizers(int p, const vector<minimizer_t> &index) 
{
	int lo = 0, hi = index.size() - 1, mid;
	while (lo <= hi) {
		mid = lo + (hi - lo) / 2;
		if (index[mid].second >= p && (!mid || index[mid - 1].second < p))
			break;
		if (index[mid].second < p) lo = mid + 1;
		else hi = mid;
	}
	assert(index[mid].second >= p);
	assert(!mid || index[mid-1].second < p);
	return mid;
}

void manual(int i, const map<int, bool> &L0, const vector<minimizer_t> &index)
{
	auto L = L0;
	int idx_i = get_minimizers(i, index), idx_j;
	int ii = idx_i;
	//print(">> {}  ", ii);
	while (index[ii].second < i+MIN_READ_SIZE) {
		//print("{:02X}  ", index[ii].first);
		bool found = (L0.find(index[ii].first) != L0.end());
		L[index[ii].first] = found;
		ii++;
	}
	auto sth = next(L.begin(), L0.size() - 1);
	int jaccard = sth->second;
	for (auto it = L.begin(); it != sth; it++)
		jaccard += it->second;
	// print (" -- {}\n", jaccard);
}

int add_sliding(auto L, auto &boundary, auto key, auto val) 
{
	int diff = 0;
	auto present = L.find(key);
	if (key <= boundary->first) {
		if (present == L.end()) {
			diff = val - boundary->second;
			boundary--;
		} else {
			diff = val - present->second;
		}
	}
	L[key] = val;
	return diff;
}

int remove_sliding(auto L, auto &boundary, auto key) 
{
	int diff = 0;
	auto present = L.find(key);
	if (present == L.end()) 
		return diff;
	if (key <= boundary->first) {
		boundary++;
		diff = int(boundary->second) - present->second;
	}
	L.erase(present);
	return diff;
}

auto refine(int x, int y, int sketch_size, const map<int, bool> &L0, const vector<minimizer_t> &index)
{
	double tau = ::tau();
	int i = x, j = x + MIN_READ_SIZE;
	print("x={}, y={}, span={}, s={}\n", x, y, y-x, sketch_size);
	
	auto L = L0;
	int idx_i = get_minimizers(i, index), idx_j;
	int ii = idx_i; //, XX=0;
	while (index[ii].second < j) {
		bool found = (L0.find(index[ii].first) != L0.end());
		L[index[ii].first] = found;
		ii++;
	}
	idx_j = ii;

	auto sth = next(L.begin(), sketch_size - 1);
	int jaccard = sth->second;
	for (auto it = L.begin(); it != sth; it++)
		jaccard += it->second;
	
	//int optI = 0, optJ = 0;
	int seed_i = -1, seed_jaccard;
	if (jaccard >= tau * sketch_size)
		seed_i = i, seed_jaccard = jaccard;
	/*if (seed_i != -1)*/ while (i < y) {
		if (index[idx_i].second < i + 1) {
			jaccard += remove_sliding(L, sth, index[idx_i].first);
			idx_i++;		
		}
		if (index[idx_j].second < j + 1) {
			bool found = L0.find(index[idx_j].first) != L0.end();
			jaccard += add_sliding(L, sth, index[idx_j].first, found);
			idx_j++;
		}

		if (jaccard >= tau * sketch_size && jaccard >= seed_jaccard) {
			seed_i = i, seed_jaccard = jaccard;
			//break;
		}
		i++, j++;
	}

/*****************************************************************************************/
	
	//candidates.push_back(make_tuple(i, j, jaccard));

 	// L has all minimizers now
 	// we are mapping p, q to i, j
 	// int p, q;
 	// int idx_p, idx_q;
 	// idx_i/idx_j points to the first/last minimizer of B_i,j
 	// idx_p/idx_q points to the first/last minimizer of A
 // 	while (true) {
 // 		if (index[idx_q].second < q + 1) {
 // 			bool found = index[idx_j].first
 // 		}
 // 		if (index[idx_j].second <= j + 1) { // should we extend?
	// 		bool found = L0.find(index[idx_j].first) != L0.end();
	// 		auto inL = L.find(index[idx_j].first);
	// 		if (index[idx_j].first <= sth->first) {
	// 			if (found && !inL->second)
	// 				jaccard += found;
	// 			if (inL == L.end()) {
	// 				jaccard -= sth->second;
	// 				sth--;
	// 			}
	// 		}
	// 		L[index[idx_j].first] = found;
	// 		idx_j++;
	// 	} 

	// 	if (jaccard >= tau * sketch_size && jaccard >= seed_jaccard) {
	// 		seed_i = i, seed_jaccard = jaccard;
	// 		break;
	// 	}
	// 	i++, j++;
	// }


	return make_pair(seed_i, seed_jaccard);
}

// s = sketch size, tau = G - delta = 1/(2exp(k*error)-1) - delta
auto search (int start, const hash_t &hash, const vector<minimizer_t> &index) 
{
	vector<int> candidates;
	map<int, bool> L0;
	uint64_t ph = 1LL << 63;

	int st = get_minimizers(start, index);
	// print(">>> {}  ", st);
	for (int mi = st; ; mi++) {
		// print("{:02X}  ", index[mi].first);

		auto &m = index[mi];
		if (m.second - start > MIN_READ_SIZE)  // < -- 1 MB?
			break; // TODO count initial kmer overlap
		if (m.first == ph) 
			continue;
		L0[m.first] = 0;
		ph = m.first;
		auto ptr = hash.find(m.first);
		if (ptr == hash.end()) continue;
		for (auto &pos: ptr->second) {
			// TODO Ignore iff hash has too many positions (high freq filter)
			// at least 1 kb spacing
			if (!(pos < start - MIN_READ_SIZE || pos > start + 2 * MIN_READ_SIZE))
				continue;
			candidates.push_back(pos);
		}
	}
	// print("\n");
	sort(candidates.begin(), candidates.end());

	int sketch_size = L0.size();
	int M = ceil(sketch_size * tau());
	vector<pair<int, int>> P;
	int px = -1, py = -1;
	for (int i = 0; i <= (int)candidates.size() - M; i++) {
		int j = i + (M - 1);
		if (candidates[j] - candidates[i] < MIN_READ_SIZE) {
			int x = max(0, candidates[j] - MIN_READ_SIZE + 1), y = candidates[i] + 1;
			if (x >= py) {
				if (px != py) P.push_back(refine(px, py, sketch_size, L0, index));
				px = x, py = y;
			} else {
				py = y;
			}
		}
	}
	if (px != py) P.push_back(refine(px, py, sketch_size, L0, index));

	return make_pair(P, double(sketch_size));
}

inline double j2md(float j, int k)
{
	double d;
	if (fabs(j) < 1e-8) d = 1;
	else if (fabs(j - 1) < 1e-8) d = 0;
	else d = (-1.0 / k) * log(2.0 * j/(1 + j));
	return 100 * (1 - d);
}


int main(void)
{
// 2	chr1	88000	121417	chr1:235525	0	+	chr1	235525	267707	32182	...	32150	31941	209	133	76	0.993499	0.992727	0.006529	0.006532	33417
// 4	chr1	91256	92392	chr1:521369	0	+	chr1	521369	522487	1118	...	1117	1092	25	18	7	0.977619	0.974130	0.022722	0.022781	1136
// 6	chr1	92387	104808	chr1:573869	0	+	chr1	573869	586415	12546	...	12359	12175	184	121	63	0.985112	0.983679	0.015038	0.015056	12421
// 10	chr1	92387	135370	chr1:224095998	0	+	chr1	224095998	224139533	43535	...	42819	42180	639	417	222	0.985077	0.983377	0.015074	0.015091	42983
// 12	chr1	92387	136258	chr1:243174377	0	+	chr1	243174377	243218157	43780	...	43446	42679	767	507	260	0.982346	0.980338	0.017865	0.017891	43871

	ifstream fin("chr1.fa");
	string l, dna;
	while (getline(fin, l)) {
		if (l[0] == '>') continue;
		dna += l;
	}
	fin.close();
	transform(dna.begin(), dna.end(), dna.begin(), ::toupper);
	dna = dna.substr(0, 4000000);

	t_start = chrono::high_resolution_clock::now();

	auto hash = get_hash(dna);
		t_end = chrono::high_resolution_clock::now();
		print("time: {:.2f}s\n", chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count() / 1000.00), t_start = t_end;

	auto mapping = search(88000, hash.second, hash.first);
		t_end = chrono::high_resolution_clock::now();
		print("time: {:.2f}s\n", chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count() / 1000.00), t_start = t_end;

	print("{}\n", mapping.first.size()); //3654?
	for (auto &p: mapping.first) /*if (p.second > .5)*/ {
		print("mapping {}: {} {:.2f}\n", p.first, p.second, j2md(p.second/mapping.second, KMER_SIZE));
	}

	return 0;
}