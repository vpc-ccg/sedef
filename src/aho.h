/// 786
/// Adapted from SCALCE source code 
/// https://github.com/sfu-compbio/scalce

/******************************************************************************/

#pragma once

/******************************************************************************/

#include <list>
#include <vector>
#include <string>
#include <cmath>

#include "common.h"
#include "search.h"

using namespace std;

/******************************************************************************/

class AHOAutomata {
	struct Trie {
		array<shared_ptr<Trie>, 4> child;
		Trie *fail, *next_to_output;
		int output;

		Trie(): output(-1), fail(0), next_to_output(0) {}

		void insert_pattern(const string &s, int level, int id);
	};

	vector<string> patterns;
	shared_ptr<Trie> trie;

private:
	void initialize_automata();

public:
	AHOAutomata();
	void search(const char *text, int len, map<int, int> &hits, int flag); // const;
};
