/// 786

/******************************************************************************/

#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <chrono>

#include <glob.h>

#include "align.h"
#include "common.h"
#include "fasta.h"
#include "hit.h"

using namespace std;

/******************************************************************************/

Hit Hit::from_bed(const string &bed)
{
	auto ss = split(bed, '\t');
	assert(ss.size() >= 10);
	
	Hit h {
		make_shared<Sequence>(ss[0], "", ss[8][0] != '+'), 0, 0, 
		make_shared<Sequence>(ss[3], "", ss[9][0] != '+'), 0, 0, 
		0, "", "", {}
	};

	h.query_start = atoi(ss[1].c_str());
	h.query_end = atoi(ss[2].c_str());

	h.ref_start = atoi(ss[4].c_str());
	h.ref_end = atoi(ss[5].c_str());

	h.name = ss[6];
	if (ss.size() >= 15) {
		h.comment = ss[14];
	}
	if (ss.size() >= 14) {
		h.jaccard = atoi(ss[13].c_str());
	}

	return h;
}

Hit Hit::from_bed(const string &bed, shared_ptr<Sequence> query, shared_ptr<Sequence> ref)
{
	auto ss = split(bed, '\t');
	assert(ss.size() >= 10);

	Hit h {
		query, 0, 0, 
		ref, 0, 0, 
		0, "", "", {}
	};

	h.query_start = atoi(ss[1].c_str());
	h.query_end = atoi(ss[2].c_str());

	h.ref_start = atoi(ss[4].c_str());
	h.ref_end = atoi(ss[5].c_str());

	assert(h.query->is_rc == (ss[8][0] != '+'));
	assert(h.ref->is_rc == (ss[9][0] != '+'));
	
	assert(!h.query->is_rc);
	if (ref->is_rc) {
		swap(h.ref_start, h.ref_end);
		h.ref_start = ref->seq.size() - h.ref_start + 1;
		h.ref_end = ref->seq.size() - h.ref_end + 1;
	}

	h.name = ss[6];
	if (ss.size() >= 15) {
		h.comment = ss[14];
	}
	if (ss.size() >= 14) {
		h.jaccard = atoi(ss[13].c_str());
	}
	if (ss.size() >= 13) {
		h.aln = Alignment::from_cigar(query->seq, ref->seq, ss[12]);
	}

	return h;
}

Hit Hit::from_wgac(const string &bed, shared_ptr<Sequence> query, shared_ptr<Sequence> ref)
{
	auto ss = split(bed, '\t');
	assert(ss.size() >= 10);
	
	Hit h {
		query, 0, 0, 
		ref, 0, 0, 
		0, "", "", {}
	};

	h.query_start = atoi(ss[1].c_str());
	h.query_end = atoi(ss[2].c_str());
	h.ref_start = atoi(ss[7].c_str());
	h.ref_end = atoi(ss[8].c_str());

	assert(h.ref->is_rc == (ss[5][0] != '+'));
	assert(!h.query->is_rc);
	if (h.ref->is_rc) {
		swap(h.ref_start, h.ref_end);
		h.ref_start = h.ref->seq.size() - h.ref_start;
		h.ref_end = h.ref->seq.size() - h.ref_end;
	}

	h.name = ss[16];
	if (ss.size() >= 15) {
		h.comment = fmt::format("match={};matchindel={}", ss[25], ss[26]);
	}

	return h;
}

/******************************************************************************/

string Hit::to_bed()
{
	assert(!query->is_rc);
	return fmt::format(
		"{}\t{}\t{}\t" // QUERY 0 1 3
		"{}\t{}\t{}\t" // REF   3 4 5
		"{}\t{}\t"     // NAME 6 SCORE 7
		"+\t{}\t"    // 8 STRAND 9
		"{}\t{}\t"     // MAXLEN 10 ALNLEN 11 
		"{}\t{}\t{}{}",  // CIGAR 12 JACCARD 13 COMMENT 14
		query->name, query_start, query_end, 
		ref->name, 
		ref->is_rc ? ref->seq.size() - ref_end + 1 : ref_start, 
		ref->is_rc ?  ref->seq.size() - ref_start + 1 : ref_end,
		name,
		aln.a.size() ? aln.error.error() : -1,
		ref->is_rc ? "-" : "+",
		// Optional fields
		// - Max. length
		// - Aln. length
		// - CIGAR
		// - Jaccard similarity
		// - Reason
		max(query_end - query_start, ref_end - ref_start),
		aln.alignment.size(),
		aln.cigar_string(),
		jaccard,
		comment,
		fmt::format(";;|>miserr={};gaperr={}", aln.error.mis_error(), aln.error.gap_error())
	);
}
