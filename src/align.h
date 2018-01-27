/// 786

/******************************************************************************/

#pragma once

/******************************************************************************/

#include <string>
#include <deque>

#include "common.h"
#include "fasta.h"
#include "hash.h"

/******************************************************************************/

typedef std::pair<int, int> Loci;
typedef list<std::pair<int, int>> Guide;

struct Hit;

struct Alignment {
	struct AlignmentError {
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

	AlignmentError error;

	void trim_front();
	void trim_back();

	static Alignment from_cigar(const std::string &a, const std::string &b, const std::string &cigar);
	static Alignment from_anchors(const std::string &qstr, const std::string &rstr,
		std::vector<Hit> &guide,
		const int side);

	std::string cigar_string();
	void cigar_from_alignment();
	void populate_nice_alignment();
	std::string print(int width=100);
	std::string print_only_alignment(int width=-1);
	void prepend_cigar(const std::deque<std::pair<char, int>> &app);
	void append_cigar(const std::deque<std::pair<char, int>> &app);

	Alignment trim();
	std::vector<Alignment> max_sum(int min_span=MIN_READ_SIZE);

	AlignmentError calculate_error();
};

/******************************************************************************/

Alignment align(const std::string &fa, const std::string &fb, 
	int match = 5, 
	int mismatch = -4, 
	int gap_open = 40, 
	int gap_extend = 1,
	int bandwidth = -1);

