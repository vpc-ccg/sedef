/// 786

/******************************************************************************/

#pragma once

/******************************************************************************/

#include <string>
#include <deque>

#include "align.h"
#include "common.h"
#include "fasta.h"
#include "hash.h"

/******************************************************************************/

struct Hit {
	/*const*/ shared_ptr<Sequence> query;
	int query_start, query_end; // query range

	/*const*/ shared_ptr<Sequence> ref;
	int ref_start, ref_end; // reference range

	int jaccard; // coordinates of seed matches
	string name, comment;

	Alignment aln;

	static Hit from_bed(const std::string &bed);
	static Hit from_bed(const std::string &bed, shared_ptr<Sequence> query, shared_ptr<Sequence> ref);
	static Hit from_wgac(const string &bed);

	std::string to_bed(bool do_rc=true);
};
