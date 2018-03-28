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

struct Hit;

struct Anchor {
	int q, r, l;
	int has_u;
};


/******************************************************************************/

struct Alignment {
	static const int MATCH = 5;
	static const int MISMATCH = -4;
	static const int GAP_OPEN = -40;
	static const int GAP_EXTEND = -1;

private:
	std::string chr_a; int start_a, end_a;
	std::string chr_b; int start_b, end_b;

	std::string a, b;
	std::string align_a, align_b, alignment;
	
	std::deque<std::pair<char, int>> cigar;

	struct AlignmentError {
		int gaps, gap_bases, mismatches, matches;	
	} error;

public:

	friend std::vector<Hit> split_alignment(Hit h);
	friend std::vector<Hit> gap_split(Hit h);
	friend bool subhit(const Hit& h, int start, int end, Hit &ho);
	friend void process(Hit hs, std::string cigar, FastaReference &reference);
	friend void trimlower(Hit &h);


	Alignment();
	Alignment(const std::string &fa, const std::string &fb);
	Alignment(const std::string &fa, const std::string &fb, const std::string &cigar);
	Alignment(const std::string &qstr, const std::string &rstr, const std::vector<Hit> &guide, const int side);
	Alignment(const string &qstr, const string &rstr, 
		const vector<Anchor> &guide, const vector<int> &guide_idx);


private: // Internal functions (modify)
	void populate_nice_alignment();
	void trim();

	void trim_front();
	void trim_back();
	
	void prepend_cigar(const std::deque<std::pair<char, int>> &app);
	void append_cigar(const std::deque<std::pair<char, int>> &app);
	void cigar_from_alignment();

	void swap();
	
public: // External functions (modify)
	void merge(Alignment &cur, const string &qstr, const string &rstr);
	// std::vector<Alignment> max_sum(int min_span=MIN_READ_SIZE);

public: // Utilities
	std::string cigar_string() const;
	std::string print(int width=100, bool only_alignment=false) const;

public:
	int span() const {
		return alignment.size();
	}
	int matches() const {
		return error.matches;
	}
	int mismatches() const {
		return error.mismatches;
	}
	int gap_bases() const {
		return error.gap_bases;
	}
	int gaps() const {
		return error.gaps;
	}
	double gap_error() const {
		return pct(error.gap_bases, error.matches + error.gap_bases + error.mismatches);
	}
	double mismatch_error() const {
		return pct(error.mismatches, error.matches + error.gap_bases + error.mismatches);
	}
	double total_error() const {
		return mismatch_error() + gap_error();
	}

public:
	// friend std::vector<Hit> generate_anchors(const std::string &query, const std::string &ref, const int kmer_size);
	friend void test(int, char** argv);
	friend void update_from_alignment(Hit &h);
	friend void stats(const std::string &ref_path, const std::string &bed_path);
};

