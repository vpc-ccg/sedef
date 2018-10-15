/// 786 

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Author: inumanag

/******************************************************************************/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <thread>
#include <mutex>

// #include <boost/asio/io_service.hpp>
// #include <boost/bind.hpp>
// #include <boost/thread/thread.hpp>

#include "common.h"
#include "search.h"
#include "search_main.h"
#include "extern/argh.h"

using namespace std;

/******************************************************************************/

extern int64_t TOTAL_ATTEMPTED;
extern int64_t JACCARD_FAILED;
extern int64_t QGRAM_NORMAL_FAILED;
extern int64_t OTHER_FAILED;
extern int64_t INTERVAL_FAILED;

/******************************************************************************/

template<typename T>
int initial_search(shared_ptr<Index> query_hash, shared_ptr<Index> ref_hash, bool is_same_genome, 
	T print_function, bool show_progress=true)
{
	Tree tree;
	int total = 0, track = 0;
	int next_to_attain = 0;

	const int TRACK_PROGRESS = 10000;
	for (int qi = 0; qi < query_hash->minimizers.size(); qi++) {
		auto &qm = query_hash->minimizers[qi];

		if (show_progress && qm.loc / TRACK_PROGRESS != track) {
			eprnn("\r |>{}<| {:.1f}% (loci={:n} hits={:n})", 
				string(int(pct(qm.loc, query_hash->seq->seq.size()) / 2) + 1, '-'), 
				pct(qm.loc, query_hash->seq->seq.size()), qm.loc, total
			);
			track = qm.loc / TRACK_PROGRESS;
		}

		if (qm.loc < next_to_attain)
			continue;

		if (Globals::Internal::DoUppercaseSeeds && qm.hash.status != Hash::Status::HAS_UPPERCASE)
			continue;

		auto hits = search(qi, query_hash, ref_hash, tree, is_same_genome,
			Globals::Search::MIN_READ_SIZE, true, false);
		int min_len = query_hash->seq->seq.size();
		for (auto &pp: hits) {
			min_len = min(min_len, pp.query_end - pp.query_start);
			print_function(pp);
		}
		total += hits.size();

		next_to_attain = (min_len >= Globals::Search::MIN_READ_SIZE 
			? qm.loc + (Globals::Search::MIN_READ_SIZE * Globals::Search::MAX_ERROR) / 2 : qm.loc);
	}
	return total;
}

/******************************************************************************/

void search_single(const string &ref_path, const string &query_chr, const string &ref_chr, bool is_ref_complement, int kmer_size, int window_size)
{
	eprn("        Parameters: READ_SIZE      = {}\n"
	     "                    MAX_ERROR      = {:.2f} ({:.2f} EDIT + {:.2f} GAP; GAPFREQ={:.3f})",
		Globals::Search::MIN_READ_SIZE, 
		Globals::Search::MAX_ERROR, 
		Globals::Search::MAX_EDIT_ERROR, 
		Globals::Search::MAX_ERROR - Globals::Search::MAX_EDIT_ERROR, 
		Globals::Search::GAP_FREQUENCY);

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
	using namespace Globals;
	argh::parser cmdl;
	cmdl.add_params({
		"-k", "--kmer", "-w", "--window", "-u", "--uppercase",
		"-e", "--error", "-E", "--edit-error", "-g", "--gap-freq", "-l", "--min-read-size",
		// "-t", "--threads", "-m", "--max-parallel-size"
	});
	cmdl.parse(argc, argv);

	cmdl({"-k", "--kmer"}, Search::KMER_SIZE) >> Search::KMER_SIZE;
	cmdl({"-w", "--window"}, Search::WINDOW_SIZE) >> Search::WINDOW_SIZE;
	cmdl({"-u", "--uppercase"}, Search::MIN_UPPERCASE) >> Search::MIN_UPPERCASE;
	cmdl({"-e", "--error"}, Search::MAX_ERROR) >> Search::MAX_ERROR;
	cmdl({"-E", "--edit-error"}, Search::MAX_EDIT_ERROR) >> Search::MAX_EDIT_ERROR;
	cmdl({"-g", "--gap-freq"}, Search::GAP_FREQUENCY) >> Search::GAP_FREQUENCY;

	// int threads = -1;
	// cmdl({"-t", "--threads"}, -1) >> threads;
	// size_t max_size = 0;
	// cmdl({"-m", "--max-parallel-size"}, 0) >> max_size;
	
	Search::MIN_READ_SIZE = KB * (1 - Search::MAX_ERROR); // 700 by default

	bool is_complement = cmdl[{"-r", "reverse"}];
	if (!cmdl(2)) {
		throw fmt::format("Not enough arguments to search");
	}

	search_single(cmdl[0], cmdl[1], cmdl[2], is_complement, 
		Search::KMER_SIZE, Search::WINDOW_SIZE);
}


