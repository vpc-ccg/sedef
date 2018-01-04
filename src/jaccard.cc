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
	if (!fin.is_open()) {
		throw fmt::format("FASTA file {} does not exist", ref_path);
	}
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

void search_main(string ref_path, string query_path, bool is_complement)
{
	auto reference = read_fasta(ref_path, is_complement);
	Index ref_hash(reference.first, reference.second);
	pair<string, string> query(reference.first, "");
	Index *query_hash = &ref_hash; 

	bool is_same_genome = (ref_path == query_path) && !is_complement;
	if (!is_same_genome) {
		query = read_fasta(query_path);
		query_hash = new Index(query.first, query.second);
	}

	eprn("Same genome:        {}", is_same_genome);
	eprn("Reverse complement: {}", is_complement);

	Tree tree;
	int total = 0, track = 0;
	for (int qi = 0; qi < query_hash->minimizers.size(); qi++) {
		auto &qm = query_hash->minimizers[qi];
		if (qm.hash.status != Hash::Status::HAS_UPPERCASE) 
			continue; // ignore N or lowercase hashes

		if (qm.loc / 10000 != track) {
			eprnn("\r   {} {:.1f}% (loci={:n} hits={:n})", 
				string(int(pct(qm.loc, query_hash->seq.size()) / 2) + 1, '-'), 
				pct(qm.loc, query_hash->seq.size()), qm.loc, total
			);
			track = qm.loc / 10000;
		}

		auto hits = search(qi, ref_hash, *query_hash, tree, is_same_genome);
		for (auto &pp: hits) {
			prn("{}", print_mapping(pp, is_complement, *query_hash, ref_hash));
		}
		total += hits.size();
	}

	eprn("");
	eprn("Total:           = {:10n}", total);
	eprn("Fails: attempts  = {:10n}\n"
	     "       Jaccard   = {:10n}\n"
	     "       interval  = {:10n}\n"
	     "       lowercase = {:10n}\n"
	     "       q-grams   = {:10n}\n"
	     "       cores     = {:10n}",
	     TOTAL_ATTEMPTED, 
	     JACCARD_FAILED, 
	     INTERVAL_FAILED, 
	     OTHER_FAILED,
	     QGRAM_NORMAL_FAILED, 
	     CORE_FAILED
	);

	// if (query_hash != &ref_hash) {
	// 	delete query_hash;
	// }
}

string print_mapping(Hit &pp, bool is_complement, 
	const Index &query, const Index &ref)
{
	if (is_complement) {
		pp.ref_start = ref.seq.size() - pp.ref_start + 1;
		pp.ref_end = ref.seq.size() - pp.ref_end + 1;
		swap(pp.ref_start, pp.ref_end);
	}
	return fmt::format(
		"{}\t{}\t{}\t{}\t{}\t{}\t\t+\t{}\t{}\t{}\t{}",
		//"{}\t{:n}\t{:n}\t{}\t{:n}\t{:n}\t\t+\t{}\t{:n}\t{}\t{}",
		query.name, pp.query_start, pp.query_end, 
		ref.name, pp.ref_start, pp.ref_end,
		is_complement ? "-" : "+",
		// Optional fields
		max(pp.query_end - pp.query_start, pp.ref_end - pp.ref_start),
		pp.jaccard,
		pp.reason
	);
}