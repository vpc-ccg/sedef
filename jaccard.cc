#include <bits/stdc++.h>
// #define FMT_HEADER_ONLY
#include <fmt/format.h>
#include <boost/math/distributions/binomial.hpp>
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

const double ERROR_RATE = 0.25;
const int WINDOW_SIZE = 16; // <-- Needs to be changed
const int MIN_READ_SIZE = 1000;

chrono::time_point<chrono::high_resolution_clock> t_start, t_end;
string reference;
unordered_map<int, int> ref_index;
unordered_map<int, int> rev_ref_index;

inline double tau(double error=ERROR_RATE)
{
	return 1 / (2 * exp(KMER_SIZE * error) - 1);
}

template<class X, class Y>
inline bool in_map(const map<X, Y> &m, X k)
{
	return m.find(k) != m.end();
}

inline double j2md(float j)
{
	double d;
	if (fabs(j) < 1e-8) d = 1;
	else if (fabs(j - 1) < 1e-8) d = 0;
	else d = (-1.0 / KMER_SIZE) * log(2.0 * j/(1 + j));
	return 100 * (1 - d);
}

inline int estM(int s, double ci=0.75)
{
	int i = ceil(s * ::tau());
	for (; i >= 0; i--)
	{
		int x = boost::math::quantile(boost::math::complement(boost::math::binomial(s, double(i) / s), (1.0 - ci) / 2));
		if (j2md(double(x) / s) >= ERROR_RATE)
			break;
	}
	return i ? i - 1 : 0;
}

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
unsigned int THRESHOLD = 1<<31;
auto get_hash (string s) 
{
	auto minimizers = get_minimizers(s);
	hash_t hash;
	for (auto &x: minimizers) {
		hash[x.first].push_back(x.second);
	}

	int ignore = (minimizers.size() * 0.001) / 100.0;
    
	map<int, int> hist;
	for (auto &m: hash)
		hist[m.second.size()] += 1;
	int sum = 0;
	for (auto i = hist.rbegin(); i != hist.rend(); i++) {
		sum += i->second;
		if (sum <= ignore) THRESHOLD = i->first;
		else break;
	}
	print("threshold {}\n", THRESHOLD);

	return make_pair(minimizers, hash);
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
		diff = boundary->second - present->second;
	}
	L.erase(present);
	return diff;
}

auto refine(
	int p, // start of read A
	int idx_p, int idx_q, // hash range for W(A)
	int x, int y, // window to check for the initial MIN_READ_SIZE match
	map<uint32_t, char> L0, // W(A)
	const vector<minimizer_t> &index)
{
	double tau = ::tau();

	int i = x, j = x + MIN_READ_SIZE;
	print("x={}, y={}, span={}, s={}\n", x, y, y-x, L0.size());
	
	auto L = L0;
	int idx_i = get_minimizers(i, index);
	int idx_j = idx_i;
	while (index[idx_j].second < j) {
		L[index[idx_j].first] = in_map(L0, index[idx_j].first), idx_j++;
	}

	auto boundary = next(L.begin(), L0.size() - 1);
	int jaccard = boundary->second;
	for (auto it = L.begin(); it != boundary; it++)
		jaccard += it->second;
	
	int seed_i = -1, seed_jaccard;
	if (jaccard >= tau * L0.size())
		seed_i = i, seed_jaccard = jaccard;
	if (seed_i != -1) while (i < y) {
		if (index[idx_i].second < i + 1) {
			jaccard += remove_sliding(L, boundary, index[idx_i].first);
			idx_i++;		
		}
		if (index[idx_j].second < j + 1) {
			jaccard += add_sliding(L, boundary, index[idx_j].first, in_map(L0, index[idx_j].first));
			idx_j++;
		}
		if (jaccard >= tau * L0.size() && jaccard >= seed_jaccard) {
			seed_i = i, seed_jaccard = jaccard;
			break;
		}
		i++, j++;
	}

/*****************************************************************************************/

 	// L has all minimizers now
 	// we are mapping p, q to i, j
 	int q = p + MIN_READ_SIZE;
 	// int idx_p, idx_q;
 	// idx_i/idx_j points to the first/last minimizer of B_i,j
 	// idx_p/idx_q points to the first/last minimizer of A
 	int prev_jaccard = jaccard;
 	vector<tuple<int, int, int, int, int, int>> candidates;
 	int last_n1 = 0, last_n2 = 0;
 	while (true) {
 		// print(">> p,q={},{} ({}) i,j={},{} ({}) jaccard={} s={} score={}\n",
 		// 	p,q,q-p,i,j,j-i,jaccard,L0.size(), j2md(double(jaccard)/L0.size(), KMER_SIZE));
 	// disallow overlaps
 		if ((i < q && p < j) || (q - p > 1000000)) {
			candidates.push_back(make_tuple(p, q - 1, i, j - 1,  j2md(double(prev_jaccard)/L0.size()), 0));
 			break;
 		}

 		if (reference[q - 1] == 'N') last_n1++; else last_n1 = 0;
 		if (reference[j - 1] == 'N') last_n2++; else last_n2 = 0;
 		if (last_n1 > 10 || last_n2 > 10) {
 			candidates.push_back(make_tuple(p, q - 1, i, j - 1,  j2md(double(prev_jaccard)/L0.size()), 1));
 			break;
 		}

 	// extend right
 		// extend A
 		if (index[idx_q].second < q + 1) {
 			auto h = index[idx_q].first;
 			if (!in_map(L0, h)) {
 				L0[h] = 0;
 				jaccard += add_sliding(L, boundary, h, in_map(L, h));
 				boundary++; // |W(A)| increases
 				jaccard += boundary->second;
 			}
 			idx_q++;
 		}
 		q++;
 		// extend Bi
 		if (index[idx_j].second < j + 1) { 
 			auto h = index[idx_j].first;
			jaccard += add_sliding(L, boundary, h, in_map(L0, h));
			idx_j++;
		} 
		j++;
		if (jaccard >= tau * L0.size()) {
			prev_jaccard = jaccard;
			continue;
		}

	// If not, remember previous and try shring from the left
		candidates.push_back(make_tuple(p, q - 1, i, j - 1,  j2md(double(prev_jaccard)/L0.size()), 2));
		// shrink A
		if (index[idx_p].second < p + 1) {
 			auto h = index[idx_p].first;
 			L0.erase(h);
 			jaccard += remove_sliding(L, boundary, h);
 			jaccard -= boundary->second;
 			boundary--;
 			idx_p++;
 		}
 		p++;
 		// shrink Bi
 		if (index[idx_i].second < i + 1) {
 			auto h = index[idx_i].first;
 			jaccard += remove_sliding(L, boundary, h);
 			idx_i++;
 		}
 		i++;
		if (jaccard >= tau * L0.size()) {
			prev_jaccard = jaccard;
			continue;
		}

	// If not, return
		break;
	}
	return candidates;
}

// s = sketch size, tau = G - delta = 1/(2exp(k*error)-1) - delta
auto search (int start, const hash_t &hash, const vector<minimizer_t> &index) 
{
	vector<int> candidates;
	map<uint32_t, char> L0;
	uint64_t ph = 1LL << 63;
	int st = get_minimizers(start, index), mi;
	for (mi = st; index[mi].second - start <= MIN_READ_SIZE; mi++) { // TODO count initial kmer overlap
		auto &m = index[mi];
		if (m.first == ph) 
			continue;
		L0[m.first] = 0;
		ph = m.first;
		auto ptr = hash.find(m.first);
		if (ptr == hash.end() || ptr->second.size() >= THRESHOLD) 
			continue;
		for (auto &pos: ptr->second) {
			// TODO Ignore iff hash has too many positions (high freq filter)
			// at least 1 kb spacing
			if (!(pos < start - MIN_READ_SIZE || pos > start + 2 * MIN_READ_SIZE))
				continue;
			candidates.push_back(pos);
		}
	}
	sort(candidates.begin(), candidates.end());
	int M = ceil(L0.size() * tau());
	print("{} candidates to look for (M={})\n", candidates.size(), M);

	vector<vector<tuple<int, int, int, int, int, int>>> P;
	int px = -1, py = -1;
	for (int i = 0; i <= (int)candidates.size() - M; i++) {
		int j = i + (M - 1);
		if (candidates[j] - candidates[i] < MIN_READ_SIZE) {
			int x = max(0, candidates[j] - MIN_READ_SIZE + 1), y = candidates[i] + 1;
			if (x >= py) {
				if (px != py) P.push_back(refine(start, st, mi, px, py, L0, index));
				px = x, py = y;
			} else {
				py = y;
			}
		}
	}
	if (px != py) P.push_back(refine(start, st, mi, px, py, L0, index));

	return P;
}

int main(void)
{
// 2	chr1	88000	121417	chr1:235525	0	+	chr1	235525	267707	32182	...	32150	31941	209	133	76	0.993499	0.992727	0.006529	0.006532	33417
// 4	chr1	91256	92392	chr1:521369	0	+	chr1	521369	522487	1118	...	1117	1092	25	18	7	0.977619	0.974130	0.022722	0.022781	1136
// *rea2 6	chr1	92387	104808	chr1:573869	0	+	chr1	573869	586415	12546	...	12359	12175	184	121	63	0.985112	0.983679	0.015038	0.015056	12421
// *rea2 10	chr1	92387	135370	chr1:224095998	0	+	chr1	224095998	224139533	43535	...	42819	42180	639	417	222	0.985077	0.983377	0.015074	0.015091	42983
// *rea2 12	chr1	92387	136258	chr1:243174377	0	+	chr1	243174377	243218157	43780	...	43446	42679	767	507	260	0.982346	0.980338	0.017865	0.017891	43871

	ifstream fin("chr1.fa");
	string l;
	string dna;
	while (getline(fin, l)) {
		if (l[0] == '>') continue;
		dna += l;
	}
	fin.close();

	//dna = dna.substr(0, 4000000);
	// for (int i = 0; i < dna.size(); i++) {
	// 	rev_ref_index[i] = reference.size();
	// 	if (isupper(dna[i])) {
	// 		ref_index[reference.size()] = i;
	// 		reference += dna[i];
	// 	}
	// }
	transform(dna.begin(), dna.end(), dna.begin(), ::toupper);
	reference = dna;

	int X,Y;
	X=88000,Y=121417; print("{}-{} ({}) ~ ", rev_ref_index[X],rev_ref_index[Y],rev_ref_index[Y]-rev_ref_index[X]);
	X=235525,Y=267707; print("{}-{} ({}) \n", rev_ref_index[X],rev_ref_index[Y],rev_ref_index[Y]-rev_ref_index[X]);
	X=91256,Y=92392; print("{}-{} ({}) ~ ", rev_ref_index[X],rev_ref_index[Y],rev_ref_index[Y]-rev_ref_index[X]);
	X=521369,Y=522487; print("{}-{} ({}) \n", rev_ref_index[X],rev_ref_index[Y],rev_ref_index[Y]-rev_ref_index[X]);

	t_start = chrono::high_resolution_clock::now();

	auto hash = get_hash(dna);
		t_end = chrono::high_resolution_clock::now();
		print("time: {:.2f}s\n", chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count() / 1000.00), t_start = t_end;

	auto pos = {88000, 91256, 92387};
	for (auto &p: pos) {
		auto mapping = search(p, hash.second, hash.first);
			t_end = chrono::high_resolution_clock::now();
			print("time: {:.2f}s\n", chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count() / 1000.00), t_start = t_end;

		print("----------------------------------\n");
		print("{} mappings in total\n", mapping.size()); //3654?
		int MA=0;
		for (auto &p: mapping) /*if (p.second > .5)*/ {
			for (auto &pp: p) {
				print("mapping: {}-{} ({}) to {}-{} ({}), id={}%, reason={}\n",
					get<0>(pp), get<1>(pp), get<1>(pp)-get<0>(pp),
					get<2>(pp), get<3>(pp), get<3>(pp)-get<2>(pp),
					get<4>(pp), get<5>(pp)
				);
			
				// auto fo = fopen(format("_S_MA{}_A.fa", MA).c_str(), "w");
				// fputs(">SEQ_A\n", fo);
				// fputs(reference.substr(get<0>(pp), get<1>(pp)-get<0>(pp)).c_str(), fo);
				// fclose(fo);
				// fo = fopen(format("_S_MA{}_B.fa", MA).c_str(), "w");
				// fputs(">SEQ_B\n", fo);
				// fputs(reference.substr(get<2>(pp), get<3>(pp)-get<2>(pp)).c_str(), fo);
				// fclose(fo);
				// auto _ = system(format("blat _S_MA{0}_A.fa _S_MA{0}_B.fa _S_out{0}.psl", MA).c_str());
				// _ = system(format("tail -n+6 _S_out{0}.psl | cut -f1,2,10-17", MA).c_str());
				// ifstream ifs(format("_S_out{0}.psl", MA).c_str());
				// getline(ifs, l);
				// getline(ifs, l);
				// getline(ifs, l);
				// getline(ifs, l);
				// getline(ifs, l);
				// getline(ifs, l);
				// int ma, as, ae, bs, be;
				// sscanf(l.c_str(), "%d %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %d %d %*s %*s %d %d",
				// 	&ma, &as, &ae, &bs, &be);

				// print("BLAT reports match={}, {}-{} to {}-{}\n", (100.0*ma)/(get<1>(pp)-get<0>(pp)), as, ae, bs, be);
			}
			MA++;
		}
		print("\n\n");
	}

	return 0;
}