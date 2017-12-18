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
	int chr_a, start_a, end_a;
	int chr_b, start_b, end_b;

	string a, b;
	string align_a, align_b, alignment;
	
	deque<pair<char, int>> cigar;
	int score;

	std::string cigar_string();
	void cigar_from_alignment();
	void populate_nice_alignment();
	std::string print(int width=100);

	alignment_t trim();
	std::vector<alignment_t> max_sum(int min_span=1000);
};

/******************************************************************************/

std::string getfasta(FastaReference &fr, const std::string &chrom, int start, int end, bool rc);
alignment_t align(std::string fa, std::string fb, 
	int match = 5, 
	int mismatch = -4, 
	int gap_open = 40, 
	int gap_extend = 1);

void align_main(std::string ref_path, std::string bed_path, int resume_after);
