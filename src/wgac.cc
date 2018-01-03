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

const int    OFF = 0;
const double MIN_ID = .95;

/******************************************************************************/

double overlap(int sa, int ea, int sb, int eb)
{
	int o = min(ea, eb) - max(sa, sb);
	// eprn("overlap: {} {} / {} {} = {}", sa, ea, sb, eb, o);
	o = max(0, o);

	return o / double(eb - sb);
}

string print_mappings(vector<Hit> &mapping, const string &ca, const string &cb, int la, int lb, string end = "\n")
{
	string out;

	sort(mapping.begin(), mapping.end(), [&](Hit pp, Hit pq) {
		auto ta = make_pair(overlap(pp.p, pp.q, OFF, OFF + la), overlap(pp.i, pp.j, OFF, OFF + lb));
		auto tb = make_pair(overlap(pq.p, pq.q, OFF, OFF + la), overlap(pq.i, pq.j, OFF, OFF + lb));
		return ta > tb;
	});
	int cnt = count_if(mapping.begin(), mapping.end(), [](Hit pp){ return pp.reason.substr(0, 2) == "OK"; });
	int prn = 0;
	double prev_score = -1;
	for (auto &pp: mapping) { 
		if (!cnt || pp.reason.substr(0, 2) == "OK") {
			double score = overlap(pp.p, pp.q, OFF, OFF + la) + overlap(pp.i, pp.j, OFF, OFF + lb);
			if (prev_score != -1 && prev_score / score > 2) break;
			out += fmt::format("   {:.2f} {:.2f} # ", 
				overlap(pp.p, pp.q, OFF, OFF + la),
				overlap(pp.i, pp.j, OFF, OFF + lb)) + print_mapping(pp, 0, ca, cb, lb) + end; 
			prev_score = score;
		}
	}
	return out;
}

/******************************************************************************/

void align_wgac(string tab_path, string ref_path)
{
	eprnn("Loading reference... ");
	unordered_map<string, string> ref;
	ifstream ff(ref_path);
	string chr, s;
	while (getline(ff, s)) {
		if (s[0] == '>') {
			chr = s.substr(1);
			eprnn("{} ", chr);
		}
		else ref[chr] += s;
	}
	eprn("done!");

	ifstream fin(tab_path.c_str());
	getline(fin, s);
	unordered_set<string> seen;
	vector<vector<string>> lines;
	while (getline(fin, s)) {
		auto ss = split(s, '\t');

		if (ss[0][3] == 'U' || ss[0].back() == 'm') continue;
		if (ss[6][3] == 'U' || ss[6].back() == 'm') continue;

		if (seen.find(ss[16]) == seen.end()) {
			seen.insert(ss[16]);
			lines.push_back(ss);
		}
	}
	
	#pragma omp parallel for
	for (int si = 0; si < lines.size(); si++) {
		auto &ss = lines[si];

		string ca = ss[0]; int sa = atoi(ss[1].c_str()), ea = atoi(ss[2].c_str());
		string cb = ss[6]; int sb = atoi(ss[7].c_str()), eb = atoi(ss[8].c_str());
		bool rb = (ss[5][0] == '_');
		auto refa = ref[ca].substr(sa - OFF, ea - sa);
		auto refb = ref[cb].substr(sb - OFF, eb - sb);
		if (rb) refb = rc(refb);

		string out;
		
		auto aln = align(refa, refb);
		aln.chr_a   = ca;
		aln.start_a = sa;
		aln.end_a   = ea;
		aln.chr_b   = cb;
		if (!rb) {
			aln.start_b = sb;
			aln.end_b   = eb;
		} else {
			aln.start_b = -eb + 1;
			aln.end_b = -sb + 1;
		}
		auto alns = aln.trim().max_sum();
		for (auto &a: alns) {
			auto err = a.calculate_error();
			out += fmt::format(
				"{}\t"
				"{}\t{}\t{}\t"
				"{}\t{}\t{}\t"
				"{}\t{:.1f}\t+\t{}\t"
				"{}\t{}\t"
				"SUB:{};{}-{};{}-{}\t",
				si,
				aln.chr_a, a.start_a, a.end_a,
				a.chr_b, 
				a.start_b < 0 ? -a.end_b + 1 : a.start_b, 
				a.end_b < 0 ? -a.start_b + 1 : a.end_b,
				ss[16], err.error(), ss[5],
				a.alignment.size(), a.cigar_string(),
				ss[17], // original size
				a.start_a - aln.start_a, a.end_a - aln.start_a,
				a.start_b - aln.start_b, a.end_b - aln.start_b
			);
			
			auto pa = a.print_only_alignment();

			auto xa = ref[ca].substr(a.start_a, a.end_a - a.start_a);
			auto xb = ref[cb].substr(a.start_b < 0 ? -a.end_b + 1 : a.start_b, a.end_b - a.start_b);
			if (rb) xb = rc(xb);
			auto x = alignment_t::from_cigar(xa, xb, a.cigar_string());
			auto px = x.print_only_alignment();

			assert(pa == px);
			out += "\n";
		}
		#pragma omp critical
		{
			prnn("{}", out);
		}
	}
}

string pax(const string &s)
{
	vector<int> lower(s.size()/1000+1, 0);
	int lowtot = 0;
	for (int i = 0; i < s.size(); i++) {
		if (islower(s[i])) lower[i / 1000]++, lowtot++;
	}

	string r=fmt::format("{}/{}", lowtot, s.size());
	for(int i = 0; i < lower.size();i++)
		r+=fmt::format(" {}:{}", i, lower[i]);
	return r;
}

void check_wgac(string bed_path, string ref_path) 
{
	unordered_map<string, string> ref;
	FastaReference fr(ref_path);

	ifstream fin(bed_path.c_str());
	string s;
	unordered_set<string> seen;
	vector<vector<string>> lines;
	while (getline(fin, s)) {
		auto ss = split(s, '\t');
		ss.erase(ss.begin());
		if (ss[0][3] == 'U' || ss[0].back() == 'm') continue;
		if (ss[3][3] == 'U' || ss[3].back() == 'm') continue;

		if (ss[0] != "chr1" || ss[3] != "chr1") continue;

		lines.push_back(ss);

		if (ref.find(ss[0]) == ref.end()) ref[ss[0]] = fr.get_sequence(ss[0]);
		if (ref.find(ss[3]) == ref.end()) ref[ss[3]] = fr.get_sequence(ss[3]);
	}
	eprn("IN: {} lines", lines.size());
	

	#define parprnn(...) out+=fmt::format(__VA_ARGS__)

	int total = 0, pass = 0, total_fails = 0;
	#pragma omp parallel for
	for (int si = 0; si < lines.size(); si++) {
		// if (si != 459) continue;

		auto &ss = lines[si];

		string ca = ss[0]; int sa = atoi(ss[1].c_str()), ea = atoi(ss[2].c_str());
		string cb = ss[3]; int sb = atoi(ss[4].c_str()), eb = atoi(ss[5].c_str());
		bool rb = (ss[9][0] == '_');

		// if (ref.find(ca) == ref.end()) ref[ca] = fr.get_sequence(ca);
		// if (ref.find(cb) == ref.end()) ref[cb] = fr.get_sequence(cb);

		auto refa = ref[ca].substr(sa - OFF, ea - sa);
		auto refb = ref[cb].substr(sb - OFF, eb - sb);
		if (rb) refb = rc(refb);

		string out;
		int i_pass = 0, i_total_fails = 0;

		auto aln = alignment_t::from_cigar(refa, refb, ss[11]);
		aln.chr_a = ca; aln.start_a = sa; aln.end_a = ea;
		aln.chr_b = cb; 
		if (!rb) aln.start_b = sb, aln.end_b   = eb;
		else aln.start_b = -eb + 1, aln.end_b = -sb + 1;
		auto alns = vector<alignment_t>{ aln };
		for (auto &a: alns) {
			auto err = a.calculate_error();
			parprnn(
				"{}\t{}\t{}\t"
				"{}\t{}\t{}\t"  
				"{}\t{}\t{}\t"
				"{}\t{:.1f}\t+\t{}\t"
				"{}\t{}\t"
				"SUB:{};{}-{};{}-{}\t",
				si,
				pax(refa), pax(refb),
				aln.chr_a, a.start_a, a.end_a,
				a.chr_b, 
				a.start_b < 0 ? -a.end_b + 1 : a.start_b, 
				a.end_b < 0 ? -a.start_b + 1 : a.end_b,
				ss[6], err.error(), ss[9],
				a.alignment.size(), "...",

				a.alignment.size(), // original size
				a.start_a - aln.start_a, a.end_a - aln.start_a, 
				a.start_b - aln.start_b, a.end_b - aln.start_b
			);

			// we need to padd it properly to assure their equality!
			// eprn("{} {}", a.a.size(), a.b.size());
			// eprn("{} {}", refa.size(), refb.size());
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
			
			Hash ha(aa);
			Hash hb(ab);

			bool success = false;

			TREE_t tree;
			vector<Hit> mappings;

			// 1. Try normal Sedef mapping
			for (int i = OFF; i < ha.seq.size(); i += 250) {
				auto m = search(i, ha, hb, tree, true);
				mappings.insert(mappings.end(), m.begin(), m.end());
				for (auto &pp: m) { 
					if (pp.reason.substr(0, 2) == "OK" 
						&& overlap(pp.p, pp.q, OFF, OFF + aa.size()) >= MIN_ID
						&& overlap(pp.i, pp.j, OFF, OFF + ab.size()) >= MIN_ID)
					{
						success = true;
						break;
					}
				}
				if (success) break;
				// break;
			}
			if (success) {
				parprnn("EXT/OK\n");
				i_pass++;
			} else {
				// eprnn("\n" + a.print());

				parprnn("\n  EXT/FAIL\n");
				parprnn("{}", print_mappings(mappings, ca, cb, aa.size(), ab.size()));

				// 2. Try full mappings without extension
				mappings = search(OFF, ha, hb, tree, true, max(aa.size(), ab.size()), false, true);
				for (auto &pp: mappings) { 
					if (pp.reason.substr(0, 2) == "OK" 
						&& overlap(pp.p, pp.q, OFF, OFF + aa.size()) >= MIN_ID
						&& overlap(pp.i, pp.j, OFF, OFF + ab.size()) >= MIN_ID)
					{
						success = true;
						break;
					}
					if (success) break;
				}
				if (success) {
					parprnn("  FULL/OK\n");
				} else {
					parprnn("  FULL/FAIL\n");
					parprnn(a.print(a.alignment.size()));
					parprnn("{}", print_mappings(mappings, ca, cb, aa.size(), ab.size()));		
				
					i_total_fails++;	
				}
			}
		}

		#pragma omp critical
		{
			/*if (i_pass == 0)*/ prnn("{}", out);
			total += alns.size();
			total_fails += i_total_fails;
			pass += i_pass;
			// if (i_pass == 0) exit(0);
		}
	}

	eprn("total:     {:>6n}\n"
		  "pass:      {:>6n} ({:.2f})\n"
		  "fail/tot:  {:>6n} ({:.2f})\n"
		  "fail/ext:  {:>6n} ({:.2f})",
		total, 
		pass, pct(pass, total),
		total_fails, pct(total_fails, total),
		total - pass  - total_fails, pct(total - pass  - total_fails, total)
	);
}
