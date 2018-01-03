/// 786

/******************************************************************************/

#pragma once

/******************************************************************************/

#include <string>
#include <deque>

#include "common.h"
#include "fasta.h"

/******************************************************************************/

struct alignment_t {
	struct aln_error_t {
		int gaps, gap_bases, mismatches, matches;

		double gap_error() {
			return pct(gap_bases, matches + gap_bases + mismatches);
		}
		double mis_error() {
			return pct(mismatches, matches + gap_bases + mismatches);
		}
		double error() {
			return mis_error() + gap_error();
		}
		double identity() {
			return 100 - error();
		}
	};

	std::string chr_a; int start_a, end_a;
	std::string chr_b; int start_b, end_b;

	std::string a, b;
	std::string align_a, align_b, alignment;
	
	std::deque<std::pair<char, int>> cigar;

	static alignment_t from_cigar(const std::string &a, const std::string &b, const std::string &cigar);
	
	std::string cigar_string();
	void cigar_from_alignment();
	void populate_nice_alignment();
	std::string print(int width=100);
	std::string print_only_alignment(int width=-1);

	alignment_t trim();
	std::vector<alignment_t> max_sum(int min_span=MIN_READ_SIZE);

	aln_error_t calculate_error();
};

/******************************************************************************/

std::string getfasta(FastaReference &fr, const std::string &chrom, int start, int end, bool rc);

alignment_t align(const std::string &fa, const std::string &fb, 
	int match = 5, 
	int mismatch = -4, 
	int gap_open = 40, 
	int gap_extend = 1);

void align_main(const std::string &ref_path, const std::string &bed_path, int resume_after);
