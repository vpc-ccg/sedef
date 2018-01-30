/// 786 

/******************************************************************************/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <bitset>
#include <unordered_set>

#include <boost/dynamic_bitset.hpp>

#include "align.h"
#include "common.h"
#include "fasta.h"
#include "hit.h"
#include "merge.h"

using namespace std;

/******************************************************************************/

void stats(const string &ref_path, const string &bed_path) 
{
	FastaReference fr(ref_path);
	ifstream fin(bed_path.c_str());
	if (!fin.is_open()) {
		throw fmt::format("BED file {} does not exist", bed_path);
	}

	string s;
	vector<pair<Hit, string>> hits;
	while (getline(fin, s)) {
		string cigar;
		Hit h = Hit::from_bed(s, &cigar);
		assert(h.ref != nullptr);
		assert(h.query != nullptr);
		if (tie(h.query->name, h.query_start, h.query_end) > tie(h.ref->name, h.ref_start, h.ref_end)) {
			swap(h.query->name, h.ref->name);
			swap(h.query_start, h.ref_start);
			swap(h.query_end, h.ref_end);
		}
		hits.push_back({h, cigar});
	}
	sort(hits.begin(), hits.end(), [](const pair<Hit, string> &a, const pair<Hit, string> &b) {
		return 
			tie(a.first.ref->is_rc, a.first.query->name, a.first.ref->name, a.first.query_start, a.first.ref_start) <
			tie(b.first.ref->is_rc, b.first.query->name, b.first.ref->name, b.first.query_start, b.first.ref_start);
	});

	int hit_count = 0, out_count = 0;
	string prev;
	for (auto &hs: hits) {
		auto &h = hs.first;
		auto &cigar = hs.second;
		string fa = fr.get_sequence(h.query->name, h.query_start, &h.query_end);
		string fb = fr.get_sequence(h.ref->name, h.ref_start, &h.ref_end);
		assert(!h.query->is_rc); 
		if (h.query->is_rc) {
			fa = rc(fa);
		}
		if (h.ref->is_rc) {
			fb = rc(fb);
		}
		assert(cigar.size());
		h.aln = Alignment(fa, fb, cigar);
		// h.aln.trim();
		// hits.push_back(h);
	// }
	// eprn("Loaded {} hits", hits.size());
	// // hits = merge(hits, 250);
	// // eprn("After merging remaining {} hits", hits.size());
	// // for (auto &h: hits) {
	// // 	prn("{}", h.to_bed(false));
	// // }
	// // exit(0);
	// for (auto &h: hits) {
		int align_length = h.aln.span();
		int indel_a = 0;
		int indel_b = 0;
		int alignB = 0;
		int matchB = 0;
		int mismatchB = 0;
		int transitionsB = 0;
		int transversionsB = 0;

		int uppercaseA = 0;
		int uppercaseB = 0;

		for (int i = 0; i < align_length; i++) {
			char a = toupper(h.aln.align_a[i]);
			char b = toupper(h.aln.align_b[i]);
			indel_a += a == '-';
			indel_b += b == '-';
			matchB += a != '-' && a == b;
			uppercaseA += (h.aln.align_a[i] != '-' && h.aln.align_a[i] != 'N' && isupper(h.aln.align_a[i]));
			uppercaseB += (h.aln.align_b[i] != '-' && h.aln.align_b[i] != 'N' && isupper(h.aln.align_b[i]));
			if (a != '-' && b != '-') {
				alignB += 1;
				if (a != b) {
					mismatchB += 1;
					if (a == 'A' || a == 'G' ) {
						transitionsB += b == 'A' || b == 'G';
						transversionsB += !(b == 'A' || b == 'G');
					} else {
						transitionsB += b == 'C' || b == 'T';
						transversionsB += !(b == 'C' || b == 'T');
					}
				}
			}
		}

		double fracMatch = double(matchB) / (alignB);
		double fracMatchIndel = double(matchB) / (align_length);
	
		double jcp = double(mismatchB) / (alignB);
		double jcK = -0.75 * log(1.0 - 4.0 / 3 * jcp);
		
		double p = double(transitionsB) / (alignB);
		double q = double(transversionsB) / (alignB);
		double w1 = 1.0 / (1 - 2.0 * p - q);
		double w2 = 1.0 / (1 - 2.0 * q);
		double k2K = 0.5 * log(w1) + 0.25 * log(w2);
		
		// TODO handle this smarter
		bool same_chr = h.query->name == h.ref->name && h.query->is_rc == h.ref->is_rc;
		int overlap = !same_chr ? 0 : max(0, 
			min(h.query_end, h.ref_end) - max(h.query_start, h.ref_start));
		bool too_big_overlap = 
			(h.query_end - h.query_start - overlap) < 100 ||
			(h.ref_end - h.ref_start - overlap) < 100;
		too_big_overlap &= same_chr;
		// too_big_overlap = 0;

		double errorScaled = (h.aln.gaps() + h.aln.mismatches()) / 
			double(h.aln.gaps() + h.aln.mismatches() + h.aln.matches());

		if (uppercaseA >= 100 && uppercaseB >= 100 && !too_big_overlap && errorScaled <= .5) {
			string l = h.to_bed(false);
			
			if (l != prev) {
				h.name = fmt::format("S{:05}", ++out_count);
				prn("{}\t"
					"{}\t{}\t{}\t{}\t{}\t{}\t"
					"{}\t{}\t{}\t{}\t{}\t{}\t"
					"{}\t{}\t{}\t"
					"{}\t{}\t{}\t{}", 
					h.to_bed(false), 
					align_length, indel_a, indel_b, alignB, matchB, mismatchB,
					transitionsB, transversionsB, fracMatch, fracMatchIndel, jcK, k2K,
					h.aln.gaps(),
					uppercaseA, uppercaseB,
					h.aln.matches(), h.aln.mismatches(), h.aln.gaps(), h.aln.gap_bases()
				);
				// if (out_count > 10) exit(0);
			}
			prev = l;
		}
		eprnn("\rProcessed hit {:n}", ++hit_count);
		// exit(0);
	}
	eprn("\nRead {:n} hits, wrote {:n} SDs", hit_count, out_count);
	eprn("\nDone!");
}

/******************************************************************************/

void get_differences(const string &ref_path, const string &bed_path,
	const string &wgac_path)
{
	map<string, boost::dynamic_bitset<>> sedef;
	map<string, boost::dynamic_bitset<>> wgac;

	string s;
	ifstream fin(bed_path);
	while (getline(fin, s)) {
		string cigar;
		Hit h = Hit::from_bed(s, &cigar);

		auto c1 = fmt::format("{}", h.query->name, "+-"[h.query->is_rc]);
		auto c2 = fmt::format("{}", h.ref->name, "+-"[h.ref->is_rc]);
		if (sedef.find(c1)==sedef.end()) sedef[c1]=boost::dynamic_bitset<>(250000000);
		if (sedef.find(c2)==sedef.end()) sedef[c2]=boost::dynamic_bitset<>(250000000);
		for (int i = h.query_start; i < h.query_end; i++) sedef[c1].set(i);
		for (int i = h.ref_start; i < h.ref_end; i++) sedef[c2].set(i);
	}

	eprn("SEDEF reading done!");

	ifstream fiw(wgac_path);
	getline(fiw, s);
	unordered_set<string> seen;
	while (getline(fiw, s)) {
		Hit h = Hit::from_wgac(s);
		auto c1 = fmt::format("{}", h.query->name, "+-"[h.query->is_rc]);
		auto c2 = fmt::format("{}", h.ref->name, "+-"[h.ref->is_rc]);		
		if (c1.size() > 6 || c2.size() > 6)
			continue;
		
		if (seen.find(h.name) == seen.end()) {
			seen.insert(h.name);
			if (wgac.find(c1)==wgac.end()) wgac[c1]=boost::dynamic_bitset<>(250000000);
			if (wgac.find(c2)==wgac.end()) wgac[c2]=boost::dynamic_bitset<>(250000000);
			for (int i = h.query_start; i < h.query_end; i++) wgac[c1].set(i);
			for (int i = h.ref_start; i < h.ref_end; i++) wgac[c2].set(i);
		}
	}

	eprn("WGAC reading done!");

	FastaReference fr(ref_path);

	int intersect = 0, wgac_only = 0, wgac_span = 0, sedef_only = 0, sedef_span = 0;

	int sedef_extra_upper = 0;
	int miss_upper = 0;

	for (auto &p: sedef) {
		auto &s = p.second;
		auto &w = wgac[p.first];

		auto seq = fr.get_sequence(p.first);

		for (int i = 0; i < seq.size(); i++) {
			if ((s[i] & (~w[i])) && isupper(seq[i]) && seq[i] != 'N') {
				sedef_extra_upper++;
			}
			if ((w[i] & (~s[i])) && isupper(seq[i]) && seq[i] != 'N') {
				miss_upper++;
			}
		}

		intersect += (s & w).count();
		wgac_only += (w & (~s)).count();
		sedef_only += (s & (~w)).count();
		sedef_span += s.count();
		wgac_span += w.count();
	}
	
	eprn("SEDEF: spans              {:12n}\n"
		 "       unique             {:12n}\n"
		 "       unique (uppercase) {:12n}\n"
		 "       misses             {:12n}\n"
		 "       misses (uppercase) {:12n}\n"
		 "WGAC:  spans              {:12n}\n"
		 "       intersects         {:12n}", 
		 sedef_span, sedef_only, sedef_extra_upper, wgac_only, 
		 miss_upper, wgac_span, intersect);
}

/******************************************************************************/

void stats_main(int argc, char **argv)
{
	if (argc < 3) {
		throw fmt::format("Not enough arguments to stats");
	}

	string command = argv[0];
	if (command == "generate") {
		stats(argv[1], argv[2]);
	} else if (command == "diff") {
		get_differences(argv[1], argv[2], argv[3]);
	} else {
		throw fmt::format("Unknown stats command");
	}
}
