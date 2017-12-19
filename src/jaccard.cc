/// 786 

/******************************************************************************/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "common.h"
#include "search.h"
#include "jaccard.h"

using namespace std;

/******************************************************************************/

extern int64_t TOTAL_ATTEMPTED;
extern int64_t JACCARD_FAILED;
extern int64_t QGRAM_NORMAL_FAILED;
extern int64_t OTHER_FAILED;
extern int64_t CORE_FAILED;
extern int64_t INTERVAL_FAILED;

/******************************************************************************/

void jaccard_search(string ref_path, string query_path, bool is_complement)
{
	string line, reference, query;

	ifstream fin;

	fin.open(ref_path.c_str());
	string ref_chr;
	while (getline(fin, line)) {
		if (line[0] != '>') reference += line;
		else ref_chr = line.substr(1);
	}
	fin.close();
	fin.clear();

	if (is_complement) {
		eprn("Reversing reference...");
		reverse(reference.begin(), reference.end());
	}

	// reference = reference.substr(0, 65000000);
	Hash ref_hash = Hash(reference);

	string query_chr = ref_chr;
	Hash *query_hash = &ref_hash; 
	if (query_path != ref_path || is_complement) {
		fin.open(query_path.c_str());
		while (getline(fin, line)) {
			if (line[0] != '>') query += line;
			else query_chr = line.substr(1);
		}
		fin.close();
		query_hash = new Hash(query);
	}

	bool allow_overlaps = (ref_path != query_path) || is_complement;
	eprn("Allowing overlaps: {}", allow_overlaps);
	eprn("Reverse complement: {}", is_complement);


	TREE_t tree;

	int total = 0;
	for (int i = 0; i < query_hash->seq.size(); i += 250) {
	// int W=51623132; for (int i = W; i < W+1; i += 250) {
		while (i < query_hash->seq.size() && query_hash->seq[i] == 'N') i++;
		while (i < query_hash->seq.size() && i % 250 != 0) i++;
		if (i % 5000 == 0) {
			double perc = 100.0 * i / double(query_hash->seq.size());
			eprnn("\r   {} {:.1f}% ({})", string(int(perc / 2) + 1, '-'), perc, i);
		}
		// prn("{}", i);

		vector<Hit> mapping = search(i, ref_hash, *query_hash, tree, allow_overlaps);
		for (auto &pp: mapping) {
			// BEDPE
			prn("{}", print_mapping(pp, is_complement, query_chr, ref_chr, ref_hash.seq.size()));
			total += 1;
		}
	}

	eprn("");
	eprn("Total:                   {:10n}", total);
	eprn("Fails: attempts          {:10n}\n"
		 "       Jaccard           {:10n}\n"
		 "       interval          {:10n}\n"
		 "       cores             {:10n}\n"
		 "       q-grams           {:10n}\n"
		 "       lowercase         {:10n}", 
		 TOTAL_ATTEMPTED, JACCARD_FAILED, INTERVAL_FAILED, 
		 CORE_FAILED, QGRAM_NORMAL_FAILED, OTHER_FAILED);
}

string print_mapping(Hit &pp, bool is_complement, 
	const string &query_chr, const string &ref_chr, int ref_size)
{
	if (is_complement) {
		pp.i = ref_size - pp.i + 1;
		pp.j = ref_size - pp.j + 1;
		swap(pp.i, pp.j);
	}
	return fmt::format("{}\t{:n}\t{:n}\t{}\t{:n}\t{:n}\t\t+\t{}\t{:n}\t{}\t{}\n",
		query_chr, pp.p, pp.q, 
		ref_chr, pp.i, pp.j,
		// pp.id, 
		is_complement ? "-" : "+",
		// Optional fields
		//int(pp.break_criteria),
		pp.q - pp.p,
		pp.jaccard,
		pp.reason
	);
}