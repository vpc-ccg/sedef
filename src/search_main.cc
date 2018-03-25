/// 786 

/******************************************************************************/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "common.h"
#include "search.h"
#include "search_main.h"

using namespace std;

/******************************************************************************/

extern int64_t TOTAL_ATTEMPTED;
extern int64_t JACCARD_FAILED;
extern int64_t QGRAM_NORMAL_FAILED;
extern int64_t OTHER_FAILED;
extern int64_t INTERVAL_FAILED;

/******************************************************************************/

bool do_uppercase = 1;
bool do_uppercase_seeds = 1;
bool do_qgram = 1;

/******************************************************************************/

void pp(shared_ptr<Index> query_hash)
{
	prn("::{} {} {}", query_hash->window_size, query_hash->kmer_size, query_hash->threshold);
	prn("::{} {} {}", query_hash->seq->is_rc, query_hash->seq->name, query_hash->seq->seq.substr(0,10));
	int64_t sum1=0,sum2=0,sum3=0;
	for (auto &m: query_hash->minimizers) {
		sum1+=m.hash.hash;
		sum2+=m.hash.status;
		sum3+=m.loc;
	}
	prn("::{} {} {}", sum1, sum2, sum3);
	sum1=sum2=sum3=0;
	for (auto &m: query_hash->index) {
		sum1+=m.first.hash;
		sum2+=m.first.status;
		for (auto x: m.second)
			sum3+=x;
	}
	prn("::{} {} {}", sum1, sum2, sum3);
}

template<typename T>
int initial_search(shared_ptr<Index> query_hash, shared_ptr<Index> ref_hash, bool is_same_genome, 
	T print_function, bool show_progress=true)
{
	Tree tree;
	int total = 0, track = 0;
	int next_to_attain = 0;

	// pp(query_hash);
	// pp(ref_hash);

	for (int qi = 0; qi < query_hash->minimizers.size(); qi++) {
		auto &qm = query_hash->minimizers[qi];

		if (show_progress && qm.loc / 10000 != track) {
			eprnn("\r |>{}<| {:.1f}% (loci={:n} hits={:n})", 
				string(int(pct(qm.loc, query_hash->seq->seq.size()) / 2) + 1, '-'), 
				pct(qm.loc, query_hash->seq->seq.size()), qm.loc, total
			);
			track = qm.loc / 10000;
		}

		if (qm.loc < next_to_attain)
			continue;

		if (do_uppercase_seeds  && qm.hash.status != Hash::Status::HAS_UPPERCASE)
			continue;
		if (!do_uppercase_seeds && qm.hash.status == Hash::Status::HAS_N)
			continue;


		// Tree tree;
		auto hits = search(qi, query_hash, ref_hash, tree, is_same_genome);
		// prn("::{} {} {}", qi, qm.loc, hits.size());
		// exit(0);
		int min_len = query_hash->seq->seq.size();
		for (auto &pp: hits) {
			min_len = min(min_len, pp.query_end - pp.query_start);
			print_function(pp);
		}
		total += hits.size();

		next_to_attain = (min_len >= MIN_READ_SIZE ? qm.loc + (MIN_READ_SIZE * MAX_ERROR) / 2 : qm.loc);
	}
	return total;
}

/******************************************************************************/

void search_parallel(const string &ref_path)
{
	FastaReference fr(ref_path);

	#ifdef _OPENMP
	eprn("Using {} threads\n", omp_get_max_threads());
	#endif

	auto T = cur_time();

	// Preload the reference
	vector<tuple<string, bool, string>> chromosomes;
	for (auto &fi: fr.index) for (int i = 0; i < 2; i++) {
		if (fi.first != "chr22" && fi.first != "chr10") 
			continue;
		chromosomes.push_back(make_tuple(fi.first, bool(i), fr.get_sequence(fi.first)));
	}
	eprn("Loading reference took {:.1f}s", elapsed(T)), T = cur_time();

	// Build index
	map<pair<string, bool>, shared_ptr<Index>> reference;
	#pragma omp parallel for schedule(dynamic)
	for (int ci = 0; ci < chromosomes.size(); ci++) {
		auto time = cur_time();

		string name, seq; bool is_rc;
		tie(name, is_rc, seq) = chromosomes[ci];

		auto idx = make_shared<Index>(make_shared<Sequence>(name, seq, is_rc));

		#pragma omp critical
		reference[{name, is_rc}] = idx;

		#pragma omp critical
		eprn("-- Built index for {:5}{} in {:5.1f}s", name, "+-"[is_rc], elapsed(time));
	}
	eprn("Building index took {:.1f}s\n", elapsed(T)), T = cur_time();
	
	// Set up instances
	vector<tuple<string, string, bool>> instances;
	for (auto &chrA: chromosomes) if (get<1>(chrA)) {
		for (auto &chrB: chromosomes) {
			if (get<2>(chrA).size() >= get<2>(chrB).size()) {
				instances.push_back(make_tuple(get<0>(chrA), get<0>(chrB), get<1>(chrB)));
			}
		}
	}
	eprn("Total {} instances", instances.size());

	int total = 0, completed = 0;
	#pragma omp parallel for schedule(dynamic)
	for (int ii = 0; ii < instances.size(); ii++) {
		auto time = cur_time();

		string chrq, chrr; bool is_rc;
		tie(chrq, chrr, is_rc) = instances[ii];

		auto itq = reference.find({chrq, false});
		assert(itq != reference.end());

		auto itr = reference.find({chrr, is_rc});
		assert(itr != reference.end());

		string out;
		int t = initial_search(itq->second, itr->second, (chrq == chrq) && !is_rc, [&](Hit &h) {
			out += h.to_bed() + "\n";
		}, false);

		#pragma omp critical
		{
			prnn("{}", out);
			fflush(stdout);
		}

		#pragma omp atomic
		total += t;

		#pragma omp atomic
		completed++;

		#pragma omp critical
		eprn("-- ({:3.0f}%%) Completed {:5}+ to {:5}{} in {:5.1f}s, {:n} hits found", 
			pct(completed, instances.size()),
			chrq, chrr, "+-"[is_rc], elapsed(time), t);
	}
	eprn("Searching took {:.1f}s", elapsed(T));

	eprn("");
	eprn("Total:           = {:10n}", total);
	eprn("Fails: attempts  = {:10n}\n"
	     "       Jaccard   = {:10n}\n"
	     "       interval  = {:10n}\n"
	     "       lowercase = {:10n}\n"
	     "       q-grams   = {:10n}\n",
	     TOTAL_ATTEMPTED, 
	     JACCARD_FAILED, 
	     INTERVAL_FAILED, 
	     OTHER_FAILED,
	     QGRAM_NORMAL_FAILED
	);
}

void search_single(const string &ref_path, const string &query_chr, const string &ref_chr, bool is_ref_complement, int kmer_size, int window_size)
{
	eprn("        Parameters: READ_SIZE      = {}\n"
	     "                    MAX_ERROR      = {:.2f} ({:.2f} EDIT + {:.2f} GAP; GAPFREQ={:.3f})",
		MIN_READ_SIZE, 
		MAX_ERROR, MAX_EDIT_ERROR, MAX_GAP_ERROR, GAP_FREQUENCY);

	bool is_same_genome = (ref_chr == query_chr) && !is_ref_complement;
	eprn("Same genome:        {}", is_same_genome);
	eprn("Reverse complement: {}", is_ref_complement);
	eprn("k-mer size:         {}", kmer_size);
	eprn("Window size:        {}", window_size);
	eprn("");

	auto T = cur_time();

	FastaReference fr(ref_path);

	string ref = fr.get_sequence(ref_chr);	
	auto ref_hash = make_shared<Index>(make_shared<Sequence>(ref_chr, ref, is_ref_complement), kmer_size, window_size);

	auto query_hash = ref_hash; 
	if (!is_same_genome) {
		string query = fr.get_sequence(query_chr);
		query_hash = make_shared<Index>(make_shared<Sequence>(query_chr, query), kmer_size, window_size);
	}
	eprn("Building index took {:.1f}s", elapsed(T)), T = cur_time();


	int total = initial_search(query_hash, ref_hash, is_same_genome, [](Hit &h) {
		prn("{}", h.to_bed());
	});

	eprn("");
	eprn("Total:           = {:10n}", total);
	eprn("Fails: attempts  = {:10n}\n"
	     "       Jaccard   = {:10n}\n"
	     "       interval  = {:10n}\n"
	     "       lowercase = {:10n}\n"
	     "       q-grams   = {:10n}",
	     TOTAL_ATTEMPTED, 
	     JACCARD_FAILED, 
	     INTERVAL_FAILED, 
	     OTHER_FAILED,
	     QGRAM_NORMAL_FAILED
	);

	exit(0);
}

/******************************************************************************/

void search_main(int argc, char **argv)
{
	do_uppercase = 1;
	do_uppercase_seeds = 1;
	do_qgram = 1;

	if (argc < 2) {
		throw fmt::format("Not enough arguments to search");
	}

	string command = argv[0];
	if (command == "single") {
		if (argc < 4) {
			throw fmt::format("Not enough arguments to align-bucket");
		}
		bool is_complement = (argc > 4 && (tolower(argv[4][0]) == 'y' || tolower(argv[4][0]) == '1'));
		int kmer_size = (argc > 5 ? atoi(argv[5]) : KMER_SIZE);
		int window_size = (argc > 6 ? atoi(argv[6]) : WINDOW_SIZE);
		search_single(argv[1], argv[2], argv[3], is_complement, kmer_size, window_size);
	} else if (command == "parallel") {
		search_parallel(argv[1]);
	} else {
		throw fmt::format("Unknown search command");
	}
}


