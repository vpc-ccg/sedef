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

pair<string, string> read_fasta(string ref_path, bool complement=false)
{
	string line, reference, chr;

	ifstream fin(ref_path);
	while (getline(fin, line)) {
		if (line[0] != '>') {
			reference += line;
		} else{
			chr = line.substr(1);
		}
	}
	fin.close();

	return make_pair(chr, complement ? rc(reference) : reference);
}

void jaccard_search(string ref_path, string query_path, bool is_complement)
{
	auto reference = read_fasta(ref_path, is_complement);
	Hash ref_hash(reference.second);
	pair<string, string> query(reference.first, "");
	Hash *query_hash = &ref_hash; 

	bool allow_overlaps = (ref_path != query_path) || is_complement;
	if (allow_overlaps) {
		query = read_fasta(query_path);
		query_hash = new Hash(query.second);
	}

	eprn("Allowing overlaps: {}", allow_overlaps);
	eprn("Reverse complement: {}", is_complement);

	TREE_t tree;
	int total = 0, track = 0;
	for (auto &qm: query_hash->minimizers) {
		auto pos = qm.second;
		if (qm.first.first) continue; // ignore N or lowercase hashes

		if (pos / 10000 != track) {
			eprnn("\r   {} {:.1f}% (loci={:n} hits={:n})", 
				string(int(pct(pos, query_hash->seq.size()) / 2) + 1, '-'), 
				pct(pos, query_hash->seq.size()), pos, total
			);
			track = pos / 10000;
		}

		auto mapping = search(pos, ref_hash, *query_hash, tree, allow_overlaps);
		for (auto &pp: mapping) {
			prn("{}", print_mapping(pp, is_complement, query.first, reference.first, ref_hash.seq.size()));
		}
		total += mapping.size();
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

	if (query_hash != &ref_hash) {
		delete query_hash;
	}
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