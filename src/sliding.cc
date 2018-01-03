/// 786

/******************************************************************************/

#include <map>

#include "common.h"
#include "sliding.h"
#include "fasta.h"

using namespace std;

/******************************************************************************/

SlidingMap::SlidingMap(): 
	query_size(0), time(0), intersection(0), limit(0)
{ 
	boundary = this->end();
} 

int SlidingMap::jaccard() 
{
	if (intersection >= limit) {
		return intersection;
	} else {
		return intersection - limit;
	}
}

bool SlidingMap::add(const hash_t &h, int BIT, int FULL) // bit: var with only one bit set
{
	// eprn("> h={}", h);

	auto it = this->lower_bound(h);  // first >=
	// v KILLING IT
	// for (; it != this->end() && it->first.first == h && (it->second & BIT); it++); 

	bool inserted = false;
	if (it != this->end() && it->first == h) { // } && !(it->second & BIT)) {
		if (it->second & BIT) return false; 
		it->second |= BIT;
	} else {
		it = this->insert({h, BIT}).first;
		inserted = true;
	}

	assert(it != this->end());
	// eprn("> it={}/{}", it->first.first, it->first.second);
	// eprn("> it={}/{}", boundary->first.first, boundary->first.second);
	// eprn("{}:{}")
	assert(boundary != this->end() || !query_size);
	if (query_size && it->first < boundary->first) {
		intersection += (it->second == FULL);
		if (inserted) {
			intersection -= (boundary->second == FULL);

			assert(boundary != this->begin()); // |S| >= 1!!!!!!!!
			boundary--;
		}
	}
	return true;
}

bool SlidingMap::remove(const hash_t &h, int BIT, int FULL) // bit: var with only one bit set
{	
	auto it = this->lower_bound(h); // first >=, and not = since...
	// assert(it != this->begin());
	// it--;
	// while (it != this->begin() && it->first.first == h && !(it->second & BIT)) it--; // find first recently added
	if (it->first != h || !(it->second & BIT)) return false;

	assert(boundary != this->end() || !query_size);
	if (query_size && it->first <= boundary->first) {
		intersection -= (it->second == FULL);
		if (it->second == BIT) { // remove
			boundary++;
			assert(boundary != this->end());
			intersection += (boundary->second == FULL);
		}
	}

	assert(it != this->end());
	assert(it != boundary);
	if (it->second == BIT) {
		this->erase(it);
	} else {
		it->second &= ~BIT;
	}
	return true;
}

void SlidingMap::add_to_query(const hash_t &h) 
{	
	if (!add(h, 1)) return;

	limit = relaxed_jaccard_estimate(++query_size);
	assert(boundary != this->end() || (query_size == 1));
	if (boundary == this->end()) boundary = this->begin();
	else boundary++;
	assert(boundary != this->end());
	intersection += (boundary->second == 3);
}

void SlidingMap::remove_from_query(const hash_t &h) 
{
	if (!remove(h, 1)) return;
	
	limit = relaxed_jaccard_estimate(--query_size);
	assert(boundary != this->begin() || !query_size);
	intersection -= (boundary->second == 3);
	if (boundary == this->begin()) boundary = this->end();
	else boundary--;
}


void SlidingMap::add_to_reference(const hash_t &h) 
{	
	if (h.first == 2) return;
	add(h, 2);
}

void SlidingMap::remove_from_reference(const hash_t &h) 
{
	if (h.first == 2) return;
	remove(h, 2);
}

/******************************************************************************/

struct TrueSlidingMap: public multimap<hash_t, bool> {
public:
	multiset<hash_t> query, ref;
	int limit;

public:
	TrueSlidingMap():limit(0) { } 
	void rewind() { }

	int jaccard() 
	{
		auto i = query.begin();
		auto j = ref.begin();

		int s = query.size();
		// if (s==0) { return 0; }
		limit = relaxed_jaccard_estimate(s);

		int total = 0;
		int intersection = 0;
		while (i != query.end() && j != ref.end() && total < s) {
			if (*i == *j) { // handle N-s?
				intersection++; 
				i++; j++;
			} else if (*i < *j) {
				i++;
			} else {
				j++;
			}
			total++;
		}

		if (intersection >= limit) {
			return intersection;
		} else {
			return intersection - limit;
		}
	}

	void add_to_query(const hash_t &h) 
	{
		query.insert(h);
	}

	void remove_from_query(const hash_t &h) 
	{
		auto it = query.find(h);
		if (it != query.end()) 
			query.erase(it);
	}

	void add_to_reference(const hash_t &h)
	{
		if (h.first) return;
		ref.insert(h);
	}

	void remove_from_reference(const hash_t &h)
	{
		if (h.first) return;
		auto it = ref.find(h);
		if (it != ref.end()) 
			ref.erase(it);
	}
};


// string prn_A(SlidingMap &winnow) 
// {
// 	string s1;
// 	s1 += fmt::format("l={} {} | ", winnow.limit, winnow.query_size);
// 	int i=0;
// 	for (auto e: winnow) {s1 += fmt::format("{}.{:x}_{} ", (int)e.first.first.first, e.first.first.second, (int)e.second); 
// 		if (i==winnow.query_size) s1+="| ";i++;}
// 		return s1;
// }

// string prn_B(TrueSlidingMap &winnow) 
// {
// 	string s1;
// 	s1 += fmt::format("l={} {} | ", winnow.limit, winnow.query.size());

// 	int x=0;
// 	auto i = winnow.query.begin(); // 1
// 	auto j = winnow.ref.begin(); // 2
// 	while (true) {
// 		if (i == winnow.query.end() && j == winnow.ref.end()) break;
// 		if (i != winnow.query.end() && j != winnow.ref.end() && *i == *j) { 
// 			s1 += fmt::format("{}.{:x}_3 ", (int)i->first, i->second);
// 			i++; j++;
// 		} else if (j == winnow.ref.end() || (i != winnow.query.end() && *i < *j)) {
// 			s1 += fmt::format("{}.{:x}_1 ", (int)i->first, i->second);
// 			i++;
// 		} else if (i == winnow.query.end() || j != winnow.ref.end()) {
// 			s1 += fmt::format("{}.{:x}_2 ", (int)j->first, j->second);
// 			j++;
// 		}
// 		if (x==winnow.query.size()) s1+="| "; x++;
// 	}

// 	return s1;
// }

// void dupert(SlidingMap &winnowA, TrueSlidingMap &winnowB, string fmt) 
// {
// 	auto E=winnowA.jaccard(), F=winnowB.jaccard();
// 	auto a=prn_A(winnowA); auto b=prn_B(winnowB);
// 	// eprn("{}\n{}\n{}", fmt, a, b);
// 	assert(a==b);
// 	if (E!=F) {
// 		eprn("jA {} != jB {}: {}", E, F, fmt);
// 		exit(1 + rand());
// 	}
// }

// void test_two_maps()
// {
// 	FastaReference fr("hg19.fa");

// 	string ref=fr.getSubSequence("chr1", 143698658, 6000);
// 	string query=fr.getSubSequence("chr1", 149239000, 6000);
	
// 	SlidingMap     winnowA;
// 	TrueSlidingMap winnowB;

// 	Hash ref_hash(ref);
// 	Hash query_hash(query);

// 	int query_start = 0, 
// 		query_end = 0,
// 		query_winnow_start = 0, 
// 		query_winnow_end = 0,
// 		ref_start = 0, 
// 		ref_end = 0,
// 		ref_winnow_start = 0, 
// 		ref_winnow_end = 0;

// 	auto fn_qe = [&]() {
// 		if (query_end >= query_hash.seq.size()) return 0;
// 		int i = 1;
// 		if (query_winnow_end < query_hash.minimizers.size() && query_hash.minimizers[query_winnow_end].second == query_end) {
// 			winnowA.add_to_query(query_hash.minimizers[query_winnow_end  ].first);
// 			winnowB.add_to_query(query_hash.minimizers[query_winnow_end++].first), i++;
// 			dupert(winnowA, winnowB, "fn_qe");
// 		}
// 		query_end++;
// 		return i;
// 	};
// 	auto fn_u_qe = [&](int i) {
// 		if (!i) return;
// 		if (i == 2) {
// 			winnowA.remove_from_query(query_hash.minimizers[--query_winnow_end].first);
// 			winnowB.remove_from_query(query_hash.minimizers[  query_winnow_end].first);
// 			dupert(winnowA, winnowB, "fn_qe/u");
// 		}
// 		query_end--;
// 	};

// 	auto fn_re = [&]() {
// 		if (ref_end >= ref_hash.seq.size()) return 0;
// 		int i = 1;
// 		if (ref_winnow_end < ref_hash.minimizers.size() && ref_hash.minimizers[ref_winnow_end].second == ref_end) {
// 			winnowA.add_to_reference(ref_hash.minimizers[ref_winnow_end  ].first);
// 			winnowB.add_to_reference(ref_hash.minimizers[ref_winnow_end++].first), i++;
// 			dupert(winnowA, winnowB, "fn_re");
// 		}
// 		ref_end++;
// 		return i;
// 	};
// 	auto fn_u_re = [&](int i) {
// 		if (!i) return;
// 		if (i == 2) {
// 			winnowA.remove_from_reference(ref_hash.minimizers[--ref_winnow_end].first);
// 			winnowB.remove_from_reference(ref_hash.minimizers[  ref_winnow_end].first);
// 			dupert(winnowA, winnowB, "fn_re/u");
// 		}
// 		ref_end--;
// 	};

// 	auto fn_qre = [&]() {
// 		if (query_end >= query_hash.seq.size() || ref_end >= ref_hash.seq.size()) return 0;
// 		return fn_qe() * 10 + fn_re();
// 	};
// 	auto fn_u_qre = [&](int i) {
// 		if (!i) return;
// 		fn_u_qe(i / 10), fn_u_re(i % 10);
// 	};

// 	auto fn_qs = [&]() {
// 		if (!query_start) return 0;
// 		int i = 1;
// 		if (query_winnow_start && query_hash.minimizers[query_winnow_start - 1].second == query_start - 1) {
// 			winnowA.add_to_query(query_hash.minimizers[--query_winnow_start].first), i++;
// 			winnowB.add_to_query(query_hash.minimizers[  query_winnow_start].first);
// 			dupert(winnowA, winnowB, "fn_qs");
// 		}
// 		query_start--;
// 		return i;
// 	};
// 	auto fn_u_qs = [&](int i) {
// 		if (!i) return;
// 		if (i == 2) {
// 			winnowA.remove_from_query(query_hash.minimizers[query_winnow_start  ].first);
// 			winnowB.remove_from_query(query_hash.minimizers[query_winnow_start++].first);
// 			dupert(winnowA, winnowB, "fn_qs/u");
// 		}
// 		query_start++;
// 	};

// 	auto fn_rs = [&]() {
// 		if (!ref_start) return 0;
// 		int i = 1;
// 		if (ref_winnow_start && ref_hash.minimizers[ref_winnow_start - 1].second == ref_start - 1) {
// 			winnowA.add_to_reference(ref_hash.minimizers[--ref_winnow_start].first), i++; 
// 			winnowB.add_to_reference(ref_hash.minimizers[  ref_winnow_start].first); 
// 			dupert(winnowA, winnowB, "fn_rs");
// 		}
// 		ref_start--;
// 		return i;
// 	};
// 	auto fn_u_rs = [&](int i) {
// 		if (!i) return;
// 		if (i == 2) {
// 			winnowA.remove_from_reference(ref_hash.minimizers[ref_winnow_start  ].first); 			
// 			winnowB.remove_from_reference(ref_hash.minimizers[ref_winnow_start++].first); 			
// 			dupert(winnowA, winnowB, "fn_rs/u");
// 		}
// 		ref_start++;
// 	};

// 	auto fn_qrs = [&]() {
// 		if (!query_start || !ref_start) return 0;
// 		return fn_qs() * 10 + fn_rs();
// 	};
// 	auto fn_u_qrs = [&](int i) {
// 		if (!i) return;
// 		fn_u_qs(i / 10), fn_u_rs(i % 10);
// 	};
	
// 	auto fns = vector<pair<function<int(void)>, function<void(int)>>>{
// 		make_pair(fn_qre, fn_u_qre),
// 		make_pair(fn_qe, fn_u_qe),
// 		make_pair(fn_re, fn_u_re),
// 		make_pair(fn_qrs, fn_u_qrs),
// 		make_pair(fn_qs, fn_u_qs),
// 		make_pair(fn_rs, fn_u_rs)
// 	};

// 	/*
// 	ext 1 to 1000
// 	ext 2 to 1000
// 	contract 2 to 500
// 	contract 1 to 500
// 	extend + 1000
// 	*/
// 	const int EXT = 1000;
// 	const int CON = EXT / 2;
// 	assert(query.size() == ref.size());
// 	assert(ref.size() == 6000);
// 	int ni = 0;
// 	ref_start = 3000, ref_end = 3000;
// 	query_start = 3000, query_end = 3000;
// 	while (ref_end < ref.size()) { // 4
// 		for (int i = 0; i < EXT; i++) { int e=fn_qre(); fn_u_qre(e); assert(e==fn_qre());  } // 3000..4000
// 		eprn("EXT/0    {}..{} {}..{}", ref_start, ref_end, query_start, query_end);
// 		for (int i = 0; i < CON; i++) { int e=fn_qe();  fn_u_qe(e);  assert(e==fn_qe());  	}  // 3000..4500
// 		eprn("EXT/1    {}..{} {}..{}", ref_start, ref_end, query_start, query_end);
// 		for (int i = 0; i < CON; i++) { int e=fn_re();  fn_u_re(e);  assert(e==fn_re());  	}  // 3000..4500
// 		eprn("EXT/2    {}..{} {}..{}", ref_start, ref_end, query_start, query_end);
// 		for (int i = 0; i < CON; i++) { int e=fn_qrs(); fn_u_qrs(e); assert(e==fn_qrs());  } // 2500...4500
// 		eprn("EXT/3    {}..{} {}..{}", ref_start, ref_end, query_start, query_end);
// 		for (int i = 0; i < CON; i++) { int e=fn_rs();  fn_u_rs(e);  assert(e==fn_rs());  	}  // 2000..4500
// 		eprn("EXT/4    {}..{} {}..{}", ref_start, ref_end, query_start, query_end);
// 		for (int i = 0; i < CON; i++) { int e=fn_qs();  fn_u_qs(e);  assert(e==fn_qs());  	}  // 2000..4500
// 		eprn("EXT/5    {}..{} {}..{}", ref_start, ref_end, query_start, query_end);
// 		ni++;
// 		if (ni > 6) exit(333);
// 		// eprn("{}..{} {}..{}", ref_start, ref_end, query_start, query_end);
// 	}
// 	eprn("exit");
// 	while (ref_start != ref_end) {
// 		if (ref_winnow_start < ref_hash.minimizers.size() && ref_hash.minimizers[ref_winnow_start].second == ref_start) {
// 			winnowA.remove_from_query(ref_hash.minimizers[ref_winnow_start].first);
// 			winnowB.remove_from_query(ref_hash.minimizers[ref_winnow_start].first);
// 			dupert(winnowA, winnowB, "X1");
// 			ref_winnow_start++;
// 		}
// 		ref_start++;
// 		if (query_winnow_start < query_hash.minimizers.size() && query_hash.minimizers[query_winnow_start].second == query_start) {
// 			winnowA.remove_from_query(query_hash.minimizers[query_winnow_start].first);
// 			winnowB.remove_from_query(query_hash.minimizers[query_winnow_start].first);
// 			dupert(winnowA, winnowB, "X2");
// 			query_winnow_start++;
// 		}
// 		query_start++;
// 	}
// 			eprn("EXT/6    {}..{} {}..{}", ref_start, ref_end, query_start, query_end);

// 	// assert(ni == 4);
// 	assert(ref_start == 6000);
// 	assert(query_start == 6000);
// 	eprn("woooorks!");
// 	exit(0);
// }
