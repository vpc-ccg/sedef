/// 786
/// Adapted from SCALCE source code 
/// https://github.com/sfu-compbio/scalce

/******************************************************************************/

#include <list>
#include <vector>
#include <string>
#include <cmath>

#include <zlib.h>

#include "common.h"
#include "aho.h"
#include "hit.h"

using namespace std;

/******************************************************************************/

// extern char _binary_patterns_bin_start;
// extern char _binary_patterns_bin_size;

/******************************************************************************/

void AHOAutomata::Trie::insert_pattern (const string &s, int level, int id) 
{
	if (level < s.size()) {
		char cx = hash_dna(s[level]);
		if (child[cx] == nullptr) {
			child[cx] = make_shared<Trie>();
		}
		child[cx]->insert_pattern(s, level + 1, id); 
	} else {
		output = id;
	}
}

AHOAutomata::AHOAutomata():
	trie(make_shared<Trie>())
{
	throw string("Not implemented");

	patterns.reserve(5000000);

	size_t dest_sz = 20 * 1000 * 1000;
	vector<uint8_t> destination(dest_sz);

	// int xx = compress( // from DeeZ
	// 	destination.data(), &dest_sz,
	// 	(uint8_t*)&_binary_patterns_bin_start, 
	// 	(int64_t)(&_binary_patterns_bin_size)
	// );
	// fwrite(destination.data(), 1, dest_sz, stdout);
	// exit(0);

	// TODO should be extern!
	char _binary_patterns_bin_start, _binary_patterns_bin_size;

	int c = uncompress(
		destination.data(), &dest_sz,
		(uint8_t*)&_binary_patterns_bin_start, 
		(int64_t)(&_binary_patterns_bin_size)
	);
	if (c != Z_OK) 
		throw string("zlib failed!");

	size_t pos = 0;
	while (1) {
		if (pos >= dest_sz)
			break;
		int16_t ln;
		memcpy(&ln, &destination[pos], sizeof(int16_t));
		pos += sizeof(int16_t);
		
		int32_t cnt;
		memcpy (&cnt, &destination[pos], sizeof(int32_t));
		pos += sizeof(int32_t);
		
		int sz = ln / 4 + (ln % 4 != 0);
		int64_t x;
		for(int i = 0; i < cnt; i++) {
			memcpy (&x, &destination[pos], sz);
			assert(sz <= 8);
			pos += sz;
			
			patterns.push_back("");
			for (int j = ln - 1; j >= 0; j--) {
				patterns.back() += "ACGT"[(x >> (2 * j)) & 3];
			}

			// prn("{}", patterns.back());

			if (patterns.back() == "AAAAAAAAAAAAAA") continue;
			if (patterns.back() == "CCCCCCCCCCCCCC") continue;
			if (patterns.back() == "GGGGGGGGGGGGGG") continue;
			if (patterns.back() == "TTTTTTTTTTTTTT") continue;

			trie->insert_pattern(patterns.back(), 0, patterns.size() - 1);
		}
	}
	initialize_automata();
}

void AHOAutomata::initialize_automata()
{
	queue<Trie*> q;
	trie->fail = trie.get();

	for (int i = 0; i < 4; i++) {
		if (trie->child[i] != nullptr) {
			trie->child[i]->fail = trie.get();
			q.push(trie->child[i].get());
		}
	}

	while (!q.empty()) {
		auto cur = q.front(); q.pop();
		for (int i = 0; i < 4; i++) {
			auto t = cur->child[i].get();
			if (t != nullptr) {
				auto f = cur->fail;
				while (f != trie.get() && f->child[i] == nullptr) {
					f = f->fail;
				}
				t->fail = (f->child[i] != nullptr ? f->child[i].get() : trie.get());
				q.push(t);
				t->next_to_output = (t->fail->output >= 0) ? t->fail : t->fail->next_to_output;
			}
		}
	}

	// This part necessitates shared_ptrs
	q.push(trie.get());
	while (!q.empty()) {
		auto cur = q.front(); q.pop();
		for (int i = 0; i < 4; i++) {
			if (cur->child[i] != nullptr) {
				q.push(cur->child[i].get());
			}
			auto tmp = cur;
			while (tmp != trie.get() && tmp->child[i] == nullptr) {
				tmp = tmp->fail;
			}
			cur->child[i] = (tmp->child[i] != nullptr ? tmp->child[i] : trie);
			tmp = cur;
			while (tmp != nullptr && tmp->output == -1) {
				tmp = tmp->next_to_output;
			}
			cur->next_to_output = tmp;
		}
	}
}

void AHOAutomata::search (const char *text, int len, map<int, int> &hits, int flag) //const
{
	auto cur = trie;
	for (int i = 0; i < len; i++) {
		cur = cur->child[hash_dna(text[i])];
		if (cur == nullptr) {
			cur = trie;
		}
		if (cur->next_to_output != nullptr) {
			hits[cur->next_to_output->output] |= flag;
		}
	}
}
