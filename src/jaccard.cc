/// 786 

/******************************************************************************/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>

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
		reference = rc(reference);
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
	int pxx = 0;
	for (auto &qm: query_hash->minimizers) {
		auto pos = qm.second;
		if (qm.first.first) continue; // ignore N hashes
		// if (pos >= 1000000) break;
	
		//while (i < query_hash->seq.size() && query_hash->seq[i] == 'N') i++;
		//while (i < query_hash->seq.size() && i % 250 != 0) i++;
		if (pos / 10000 != pxx) eprnn("\r   {} {:.1f}% (loci={:n} hits={:n})", 
			string(int(pct(pos, query_hash->seq.size()) / 2) + 1, '-'), 
			pct(pos, query_hash->seq.size()), pos, total), pxx = pos / 10000;

		auto mapping = search(pos, ref_hash, *query_hash, tree, allow_overlaps);
		for (auto &pp: mapping)
			prn("{}", print_mapping(pp, is_complement, query_chr, ref_chr, ref_hash.seq.size()));
		// prn("// {} ~ {} ended", qm.first, pos);
		total += mapping.size();

		// eprn("--- {}s", 
		// chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - time).count() / 1000.00);
		// time = chrono::high_resolution_clock::now() ;
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
	return fmt::format(
		"{}\t{}\t{}\t{}\t{}\t{}\t\t+\t{}\t{}\t{}\t{}",
		//"{}\t{:n}\t{:n}\t{}\t{:n}\t{:n}\t\t+\t{}\t{:n}\t{}\t{}",
		query_chr, pp.p, pp.q, 
		ref_chr, pp.i, pp.j,
		// pp.id, 
		is_complement ? "-" : "+",
		// Optional fields
		//int(pp.break_criteria),
		max(pp.q - pp.p, pp.j - pp.i),
		pp.jaccard,
		pp.reason
	);
}