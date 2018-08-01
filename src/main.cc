/// 786 

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Author: inumanag

/******************************************************************************/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <locale>

#include "search_main.h"
#include "align_main.h"
#include "stats_main.h"
#include "merge.h"

using namespace std;

/******************************************************************************/

void print_help()
{
	eprn(
		"Basic usage: sedef [command] [args...]\n"
		"All meaningful output is redirected to stdout. \n \n"
		"It's recommended to use sedef.sh unless you really know what you are doing. "
		"List of commands: \n"
		"- search [reference.fa] [chrQuery] [chrRef] \n"
		"  finds all initial SDs by aligning chrQuery onto the chrRef. \n"
                "  reference.fa should be soft-masked (i.e. all common repeats are in lowercase characters. \n"
		"  params: \n"
		"    -r, --reverse \n"
		"      set if you want reference to be reverse-complemented \n"
		"    -k, --kmer [kmer] \n"
		"      set index k-mer size (default: 12) \n"
		"    -w, --window \n"
		"      set winnowing window size (default: 16) \n"
		"    -u, --uppercase \n"
		"      set minimum number of uppercase (non-masked) letters in the seed hit (default: 12) \n"
		"    -e, --error \n"
		"      set maximum divergence rate between SDs (default: 0.3) \n"
		"    -E, --edit-error \n"
		"      set maximum small mutation rate (default: 0.15) \n"
		"    -g, --gap-freq \n"
		"      set gap frequency rate (default: 0.005) \n"
		"- align bucket -n [count] [bed_directory(/)] [buckets/] \n"
		"  bucket BEDs from [bed_directory] directory (or file) in the [count] bins (files) \n"
		"  in [buckets/] for balanced parallel alignment stage \n"
		"  params: \n"
		"    -n, --bins \n"
		"      number of bins \n"
		"- align generate [genome.fa] [initial.bed] \n"
		"  generates true alignments for [initial.bed] \n"
		"  WARNING: takes a lot of time--- run in parallel via sedef.sh if possible. \n"
		"  params: \n"
		"    -k, --kmer \n"
		"      set seeding k-mer size (not necessarily equal to search k-mer size) \n"
		"    --match, --mismatch, --gap-open, --gap-extend \n"
		"      set alignment parameters for initial seed alignment (default: 5, -4, -40, -1) \n"
		"    --extend-ratio, --max-extend \n"
		"      set parameters for seed SD heuristic extension (default: 5, 15000) \n"
		"      formula: each SD of len l is extended on both sides with \n"
		"               min(max-extend, extend-ratio * l) bases \n"
		"    --merge-dist \n"
		"      set merge distance threshold (default: 250) \n"
		"- stats generate [genome.fa] [input.bed] \n"
		"  generates and filters final SD alignments for [input.bed] \n"
		"  params: \n"
		"    --max-ok-gap \n"
		"      maximal allowed gap size expressed as percentage of SD length (default: 50) \n"
		"    --min-split \n"
		"      minimal allowed SD length after splitting large gaps (default: 1000) \n"
		"    --uppercase \n"
		"      minimum amount of uppercase (non-masked) characters in SD (default: 100) \n"
		"    --max-error \n"
		"      maximum wgac-scaled error rate allowed for SD (default: 0.5) \n"
		"- help \n"
		"  Displays this help message \n"
		"\n"
		"Questions? Please contact inumanag at mit dot edu."
	);
}

/******************************************************************************/

int main(int argc, char **argv)
{
	ios_base::sync_with_stdio(0);
	setlocale(LC_NUMERIC, "en_US.UTF-8");
	if (argc < 2) {
		eprn("Arguments missing: please run sedef help for more information.");
		exit(1);
	}

	eprnn("🍁  🐚    SEDEF {}; arguments: ", string(GITVER) == "" ? "vpc" : GITVER);
	#ifdef __SSE4_1__
		eprnn(" (SSE4.1)");
	#elif defined __SSE2__
		eprnn(" (SSE2)");
	#else
		eprnn(" (No SSE)");
	#endif
	for (int i = 0; i < argc; i++) {
		eprnn(" {}", argv[i]);
	}
	eprn("");

	string command = argv[1];

	if (argc < 3 && command != "help") {
		eprn("Arguments missing: please run sedef help for more information.");
		exit(1);
	}

	
	try {
		if (command == "help") {
			print_help();
			exit(0);
		} else if (command == "search") {
			search_main(argc - 2, argv + 2);
		} else if (command == "align") {
			align_main(argc - 2, argv + 2);
		} else if (command == "stats") {
			stats_main(argc - 2, argv + 2);
		} else {
			eprn("Whoops, invalid command!");
		}
	} catch (string &s) {
		eprn("Error: {}", s);
		eprn("Double-check the parameters: run sedef --help for explanation.");
		exit(1);
	} catch (exception &e) {
		eprn("Error: {}", e.what());
		exit(1);
	} 
	
	return 0;
}
