/// 786

/******************************************************************************/

#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <thread>
#include <unordered_set>
#include <unordered_map>

#include "common.h"
#include "search.h"
#include "jaccard.h"
#include "fasta.h"
#include "filter.h"
#include "align.h"
#include "aho.h"
#include "sliding.h"

using namespace std;

/******************************************************************************/

const double MIN_ID = .95; // Minimum overlap percentage to call successful WGAC hit

/******************************************************************************/

double overlap(int sa, int ea, int sb, int eb)
{
	int o = max(0, min(ea, eb) - max(sa, sb));
	return o / double(eb - sb);
}

string print_mappings(vector<Hit> &mapping, const Index &ca, const Index &cb, int la, int lb, string end = "\n")
{
	string out;

	sort(mapping.begin(), mapping.end(), [&](Hit pp, Hit pq) {
		auto ta = make_pair(overlap(pp.query_start, pp.query_end, 0, la), overlap(pp.ref_start, pp.ref_end, 0, lb));
		auto tb = make_pair(overlap(pq.query_start, pq.query_end, 0, la), overlap(pq.ref_start, pq.ref_end, 0, lb));
		return ta > tb;
	});
	int cnt = count_if(mapping.begin(), mapping.end(), [](Hit pp){ return pp.reason.substr(0, 2) == "OK"; });
	int prn = 0;
	double prev_score = -1;
	for (auto &pp: mapping) { 
		if (!cnt || pp.reason.substr(0, 2) == "OK") {
			double score = overlap(pp.query_start, pp.query_end, 0, la) 
				+ overlap(pp.ref_start, pp.ref_end, 0, lb);
			if (prev_score != -1 && prev_score / score > 2) break;
			out += fmt::format("   {:.2f} {:.2f} # ", 
				overlap(pp.query_start, pp.query_end, 0, la),
				overlap(pp.ref_start, pp.ref_end, 0, lb)) + print_mapping(pp, 0, ca, cb) + end; 
			prev_score = score;
		}
	}
	return out;
}

/******************************************************************************/

void align_wgac(string ref_path, string tab_path)
{
	eprnn("Loading reference... ");
	unordered_map<string, string> ref;
	FastaReference fr(ref_path);

	ifstream fin(tab_path.c_str());
	if (!fin.is_open()) {
		throw fmt::format("WGAC TAB file {} does not exist", tab_path);
	}

	string s;
	getline(fin, s); // header
	unordered_set<string> seen;
	vector<vector<string>> lines;
	while (getline(fin, s)) {
		auto ss = split(s, '\t');

		if (ss[0][3] == 'U' || ss[0].back() == 'm') 
			continue;
		if (ss[6][3] == 'U' || ss[6].back() == 'm') 
			continue;

		if (seen.find(ss[16]) == seen.end()) {
			seen.insert(ss[16]);
			lines.push_back(ss);
		}

		if (ref.find(ss[0]) == ref.end()) 
			ref[ss[0]] = fr.get_sequence(ss[0]);
		if (ref.find(ss[3]) == ref.end()) 
			ref[ss[3]] = fr.get_sequence(ss[3]);
	}
	
	prnn("chromA\tstartA\tendA\tchromB\tstartB\tendB\tname\tdiff\tstrandA\tstrandB\talnSize\tcigar\twgacAlnSize\twgacLocation");
	#pragma omp parallel for
	for (int si = 0; si < lines.size(); si++) {
		auto &ss = lines[si];

		string ca = ss[0]; int sa = atoi(ss[1].c_str()), ea = atoi(ss[2].c_str());
		string cb = ss[6]; int sb = atoi(ss[7].c_str()), eb = atoi(ss[8].c_str());
		bool rb = (ss[5][0] == '_');
		auto refa = ref[ca].substr(sa, ea - sa);
		auto refb = ref[cb].substr(sb, eb - sb);
		if (rb) refb = rc(refb);

		string out;
		
		auto aln = align(refa, refb);
		aln.chr_a = ca; aln.start_a = sa; aln.end_a = ea;		
		aln.chr_b = cb;
		if (!rb) aln.start_b = sb, aln.end_b = eb;
		else aln.start_b = -eb + 1, aln.end_b = -sb + 1;
		auto a = aln.trim();
		out += fmt::format(
			// "{}\t"
			"{}\t{}\t{}\t"
			"{}\t{}\t{}\t"
			"{}\t{:.1f}\t+\t{}\t"
			"{}\t{}\t"
			"SUB:{};{}-{};{}-{}\n",
			// si,
			aln.chr_a, a.start_a, a.end_a,
			a.chr_b, 
			a.start_b < 0 ? -a.end_b + 1 : a.start_b, 
			a.end_b < 0 ? -a.start_b + 1 : a.end_b,
			ss[16], a.calculate_error().error(), ss[5],
			a.alignment.size(), a.cigar_string(),
			ss[17], // original size
			a.start_a - aln.start_a, a.end_a - aln.start_a,
			a.start_b - aln.start_b, a.end_b - aln.start_b
		);
		
		#pragma omp critical 
		prnn("{}", out);
	}
}

void check_wgac(string ref_path, string bed_path) 
{
	unordered_map<string, string> ref;
	FastaReference fr(ref_path);

	ifstream fin(bed_path.c_str());
	if (!fin.is_open()) {
		throw fmt::format("BED file {} does not exist", bed_path);
	}

	vector<vector<string>> lines;
	string s;
	while (getline(fin, s)) {
		auto ss = split(s, '\t');
		ss.erase(ss.begin());

		if (ss[0][3] == 'U' || ss[0].back() == 'm') 
			continue;
		if (ss[3][3] == 'U' || ss[3].back() == 'm') 
			continue;
		if (ss[0] != "chr1" || ss[3] != "chr1") 
			continue;

		lines.push_back(ss);

		if (ref.find(ss[0]) == ref.end()) 
			ref[ss[0]] = fr.get_sequence(ss[0]);
		if (ref.find(ss[3]) == ref.end()) 
			ref[ss[3]] = fr.get_sequence(ss[3]);
	}
	eprn("Loaded {} lines", lines.size());

	
	#define parprnn(...) out += fmt::format(__VA_ARGS__)

	int total = 0, pass = 0, total_fails = 0;

	prnn("chromA\tstartA\tendA\tchromB\tstartB\tendB\tname\tdiff\tstrandA\tstrandB\talnSize\tcigar\twgacAlnSize\twgacLocation");
	#pragma omp parallel for
	for (int si = 0; si < lines.size(); si++) {
		auto &ss = lines[si];

		string ca = ss[0]; int sa = atoi(ss[1].c_str()), ea = atoi(ss[2].c_str());
		string cb = ss[3]; int sb = atoi(ss[4].c_str()), eb = atoi(ss[5].c_str());
		bool rb = (ss[9][0] == '_');

		auto refa = ref[ca].substr(sa, ea - sa);
		auto refb = ref[cb].substr(sb, eb - sb);
		if (rb) refb = rc(refb);

		string out;
		int i_pass = 0, i_total_fails = 0;

		auto a = alignment_t::from_cigar(refa, refb, ss[11]);
		a.chr_a = ca; a.start_a = sa; a.end_a = ea;
		a.chr_b = cb; 
		if (!rb) a.start_b = sb, a.end_b = eb;
		else a.start_b = -eb + 1, a.end_b = -sb + 1;
		
		auto err = a.calculate_error();
		parprnn(
			"{}\t{}\t{}\t"  
			"{}\t{}\t{}\t"
			"{}\t{:.1f}\t+\t{}\t"
			"{}\t{}\t"
			"{}\t",
			a.chr_a, a.start_a, a.end_a,
			a.chr_b, 
			a.start_b < 0 ? -a.end_b + 1 : a.start_b, 
			a.end_b < 0 ? -a.start_b + 1 : a.end_b,
			ss[6], err.error(), ss[9],
			a.alignment.size(),
			ss[12], ss[13]
		);

		// we need to padd it properly to assure their equality!
		string aa = a.a, ab = a.b;
		if (aa.size() < ab.size()) {
			aa += ref[ca].substr(a.end_a, ab.size() - aa.size());
		} else if (aa.size() > ab.size()) {
			if (!rb) {
				ab += ref[cb].substr(a.end_b, aa.size() - ab.size());
			} else {
				int real_start = sb + (refb.size() - (a.end_b - a.start_b));
				string r = ref[cb].substr(real_start - (aa.size() - ab.size()), aa.size() - ab.size());
				ab += rc(r);
			}
		}
		assert(aa.size() == ab.size());

		Index ha("A", aa), hb("B", ab);
		vector<Hit> mappings;
		Tree tree;
		bool success = false;

		// 1. Try normal Sedef mapping
		for (int i = 0; i < hb.minimizers.size(); i++) {
			auto &qm = hb.minimizers[i];
			if (qm.hash.status != Hash::Status::HAS_UPPERCASE) 
				continue;
			auto m = search(i, ha, hb, tree, true);
			mappings.insert(mappings.end(), m.begin(), m.end());
			for (auto &pp: m) { 
				if (pp.reason.substr(0, 2) == "OK" && 
					overlap(pp.query_start, pp.query_end, 0, aa.size()) >= MIN_ID && 
					overlap(pp.ref_start, pp.ref_end, 0, ab.size()) >= MIN_ID)
				{
					success = true;
					break;
				}
			}
			if (success) {
				parprnn("EXT/OK\n");
				i_pass++;
				break;
			}
		}
		if (!success) {
			parprnn("\n  EXT/FAIL\n");
			parprnn("{}", print_mappings(mappings, ha, hb, aa.size(), ab.size()));

			// 2. Try full mappings without extension
			mappings = search(0, ha, hb, tree, true, max(aa.size(), ab.size()), false, true);
			for (auto &pp: mappings) { 
				if (pp.reason.substr(0, 2) == "OK" 
					&& overlap(pp.query_start, pp.query_end, 0, aa.size()) >= MIN_ID
					&& overlap(pp.ref_start, pp.ref_end, 0, ab.size()) >= MIN_ID)
				{
					success = true;
					break;
				}
				if (success) {
					parprnn("  FULL/OK\n");
					break;
				}	
			}
			if (!success) {
				parprnn("  FULL/FAIL\n");
				parprnn(a.print(a.alignment.size()));
				parprnn("{}", print_mappings(mappings, ha, hb, aa.size(), ab.size()));
				i_total_fails++;	
			}
		}

		#pragma omp critical
		{
			if (i_pass == 0) prnn("{}", out);
			total++;
			total_fails += i_total_fails;
			pass += i_pass;
		}
	}

	eprn("Total:        {:>6n}\n"
	     "Pass:         {:>6n} ({:.2f})\n"
	     "Fail [extn]:  {:>6n} ({:.2f})"
	     "Fail [full]:  {:>6n} ({:.2f})\n",
		total, 
		pass, pct(pass, total),
		total_fails, pct(total_fails, total),
		total - pass  - total_fails, pct(total - pass  - total_fails, total)
	);
}

/******************************************************************************/

void wgac_main (int argc, char **argv)
{
	if (argc < 3) {
		throw fmt::format("Not enough arguments to wgac");
	}

	string command = argv[0];
	if (command == "align") {
		align_wgac(argv[1], argv[2]);
	} else if (command == "check") {
		check_wgac(argv[1], argv[2]);
	} else {
		throw fmt::format("Unknown wgac command");
	}
}

