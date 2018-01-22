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
	if (ss.size() >= 14) {
		h.comment = ss[13];
	}
	if (ss.size() >= 13) {
		h.aln = Alignment::from_cigar(query->seq, ref->seq, ss[12]);
	}

	return h;
}

Hit Hit::from_wgac(const string &bed)
{
	auto ss = split(bed, '\t');
	assert(ss.size() >= 27);
	
	Hit h {
		make_shared<Sequence>(ss[0], "", false), 
		atoi(ss[1].c_str()), 
		atoi(ss[2].c_str()),
		make_shared<Sequence>(ss[6], "", ss[5][0] != '+'), 
		atoi(ss[7].c_str()), 
		atoi(ss[8].c_str()),
		0, 
		ss[16], 
		fmt::format("err={:.1f}", 100 - 100 * atof(ss[26].c_str())), 
		{}
	};

	assert(h.ref->is_rc == (ss[5][0] != '+'));
	assert(!h.query->is_rc);

	return h;
}

/******************************************************************************/

string Hit::to_bed(bool do_rc)
{
	assert(!query->is_rc);
	return fmt::format(
		"{}\t{}\t{}\t" // QUERY 0 1 3
		"{}\t{}\t{}\t" // REF   3 4 5
		"{}\t{}\t"     // NAME 6 SCORE 7
		"{}\t{}\t"     // STRAND 8 STRAND 9
		"{}\t{}\t"     // MAXLEN 10 ALNLEN 11 
		"{}\t{}",      // CIGAR 12 COMMENT 13
		query->name, query_start, query_end, 
		ref->name, 
		do_rc && ref->is_rc ? ref->seq.size() - ref_end + 1 : ref_start, 
		do_rc && ref->is_rc ?  ref->seq.size() - ref_start + 1 : ref_end,
		name,
		aln.a.size() ? fmt::format("{:.1f}", aln.error.error()) : "",
		query->is_rc ? "-" : "+",
		ref->is_rc ? "-" : "+",
		// Optional fields
		// - Max. span
		// - Aln. length
		// - CIGAR
		// - Jaccard similarity
		// - Reason
		max(query_end - query_start, ref_end - ref_start),
		aln.alignment.size(),
		aln.cigar_string(),
		fmt::format("{}\t{}", 
			aln.a.size() ? fmt::format("m={:.1f}\tg={:.1f}", aln.error.mis_error(), aln.error.gap_error()) : "",
			comment.size() ? comment + ";" : ""
		)
	);
}
