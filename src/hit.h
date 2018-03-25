/// 786

/******************************************************************************/

#pragma once

/******************************************************************************/

#include <string>
#include <deque>

#include "align.h"
#include "common.h"
#include "fasta.h"

/******************************************************************************/

struct Hit {
	/*const*/ shared_ptr<Sequence> query;
	int query_start, query_end; // query range

	/*const*/ shared_ptr<Sequence> ref;
	int ref_start, ref_end; // reference range

	int jaccard; // coordinates of seed matches
	string name, comment;

	Alignment aln;

public:
	static Hit from_bed(const std::string &bed, std::string *cigar = nullptr);
	static Hit from_bed(const std::string &bed, shared_ptr<Sequence> query, shared_ptr<Sequence> ref);
	static Hit from_wgac(const string &bed);

public:
	std::string to_bed(bool do_rc=true, bool with_cigar=true) const;
	bool operator<(const Hit &h) const {
		return 
			tie(query_start, query_end, ref_start, ref_end) < 
			tie(h.query_start, h.query_end, h.ref_start, h.ref_end);
	}

public:
	void extend(const double factor, const int max_extend);
};

void update_from_alignment(Hit &h);
