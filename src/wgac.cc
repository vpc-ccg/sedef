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

string print_mappings(vector<Hit> &mapping, const string &ca, const string &cb, int la, int lb)
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
			out += fmt::format("<< {:.2f} {:.2f} ~ ", 
				overlap(pp.p, pp.q, OFF, OFF + la),
				overlap(pp.i, pp.j, OFF, OFF + lb)) + print_mapping(pp, 0, ca, cb, lb) + ">> "; 
			prev_score = score;
		}
	}
	return out;
}

void check_wgac(string bed_path, string ref_path) 
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

	ifstream fin(bed_path.c_str());
	getline(fin, s);
	unordered_set<string> seen;
	vector<vector<string>> lines;
	while (getline(fin, s)) {
		auto ss = split(s, '\t');

		if (ss[0][3] == 'U' || ss[0].back() == 'm') continue;
		if (ss[6][3] == 'U' || ss[6].back() == 'm') continue;
		// if (ss[0] != "chr1" || ss[6] != "chr1") continue;

		if (seen.find(ss[16]) == seen.end()) {
			seen.insert(ss[16]);
			lines.push_back(ss);
		}
	}
	// vector<vector<string>> lines = 

	// auto xx =  set<int>{980,990};
	// 	408,
	// 	543,
	// 	136,
	// 	272,
	// 	678,
	// 	813,
	// 	948
	// };
	
	int total = 0, pass = 0, total_fails = 0;
	// #pragma omp parallel for
	for (int si = 0; si < lines.size(); si++) {
		// if (xx.find(si) == xx.end()) continue;
		// #pragma omp critical
		// {
		// 	eprn("{}~~", si);
		// }

		auto &ss = lines[si];

		string ca = ss[0]; int sa = atoi(ss[1].c_str()), ea = atoi(ss[2].c_str());
		string cb = ss[6]; int sb = atoi(ss[7].c_str()), eb = atoi(ss[8].c_str());
		bool rb = (ss[5][0] == '_');

		auto refa = ref[ca].substr(sa - OFF, ea - sa);
		auto refb = ref[cb].substr(sb - OFF, eb - sb);
		if (rb) {
			reverse(refb.begin(), refb.end());
			transform(refb.begin(), refb.end(), refb.begin(), rev_dna);
		}

		string out;
		// out += fmt::format(
		// 	"{}\t{:n}\t{:n}\t{}\t{:n}\t{:n}\t"
		// 	"{}\t{}\t+\t{}" 
		// 	ss[0], atoi(ss[1].c_str()), atoi(ss[2].c_str()), 
		// 	ss[6], atoi(ss[7].c_str()), atoi(ss[8].c_str()), 
		// 	ss[16],  ss[5],
		// 	atoi(ss[17].c_str()), atof(ss[25].c_str()), atof(ss[26].c_str()), 
		// );

		int i_pass = 0, i_total_fails = 0;

		auto aln = align(refa, refb);
		aln.chr_a   = ca;
		aln.start_a = sa;
		aln.end_a   = ea;
		aln.chr_b   = cb;
		aln.start_b = sb;
		aln.end_b   = eb;
		auto alns = aln.trim().max_sum();
		for (auto &a: alns) {
			auto err = a.calculate_error();
			out += fmt::format(
				"{}\t{}\t{}\t"  
				"{}\t{}\t{}\t"
				"{}\t{:.1f}\t+\t{}\t"
				"{}\t{}\t"
				"SUB:{};{}-{};{}-{}\t",
				aln.chr_a, a.start_a, a.end_a,
				a.chr_b, a.start_b, a.end_b, 
				ss[16], err.error(), ss[5],
				a.alignment.size(), a.cigar_string(),

				ss[17], // original size
				a.start_a - aln.start_a, a.end_a - aln.start_a, 
				a.start_b - aln.start_b, a.end_b - aln.start_b
			);

			// we need to padd it properly to assure their equality!
			// eprn("{} {}", a.a.size(), a.b.size());
			// eprn("{} {}", refa.size(), refb.size());
			string aa = a.a, ab = a.b;
			if (aa.size() < ab.size()) {
				aa += ref[ca].substr(a.end_a, ab.size() - aa.size());
			}
			if (aa.size() > ab.size()) {
				if (!rb) {
					ab += ref[cb].substr(a.end_b, aa.size() - ab.size());
				} else {
					int real_start = sb + (refb.size() - (a.end_b - a.start_b));
					string rc = ref[cb].substr(real_start - (aa.size() - ab.size()), aa.size() - ab.size());
					reverse(rc.begin(), rc.end());
					transform(rc.begin(), rc.end(), rc.begin(), rev_dna);
					ab += rc;
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
				auto m = search(i, ha, hb, tree, true, MIN_READ_SIZE);
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
			}
			if (success) {
				out += fmt::format("EXT/OK\n");
				i_pass++;
				continue;
			}
			out += fmt::format("EXT/FAIL;");
			// out += print_mappings(mappings, ca, cb, a.a.size(), a.b.size());

			// 2. Try full mappings without extension
			mappings = search(OFF, ha, hb, tree, true, max(aa.size(), ab.size()));
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
				out += fmt::format("FULL/OK\n");
				continue;
			}
			out += fmt::format("FULL/FAIL\n");
			//out += a.print();
			//out += print_mappings(mappings, ca, cb, a.a.size(), a.b.size());		
			
			i_total_fails++;	
		}

		#pragma omp critical
		{
			prnn("{}", out);
			total += alns.size();
			total_fails += i_total_fails;
			pass += i_pass;
		}
		// if (i_pass != alns.size()) 			exit(0);
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
