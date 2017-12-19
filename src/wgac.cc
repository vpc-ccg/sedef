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
			out += "\t\t" + fmt::format("{:.2f} {:.2f} ~ ", 
				overlap(pp.p, pp.q, OFF, OFF + la),
				overlap(pp.i, pp.j, OFF, OFF + lb)) + print_mapping(pp, 0, ca, cb, lb); 
			prev_score = score;
		}
	}
	return out;
}

void check_wgac(string bed_path) 
{
	eprnn("Loading reference... ");
	unordered_map<string, string> ref;
	ifstream ff("data/hg19/chr1.fa");
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
		if (ss[0] != "chr1" || ss[6] != "chr1") continue;
		if (seen.find(ss[16]) == seen.end()) {
			seen.insert(ss[16]);
			lines.push_back(ss);
		}
	}
	// vector<vector<string>> lines = 

	// auto xx =  set<int>{0,
	// 	408,
	// 	543,
	// 	136,
	// 	272,
	// 	678,
	// 	813,
	// 	948
	// };
	
	int total = 0, pass = 0, total_fails = 0;
	#pragma omp parallel for
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
		out += fmt::format("{:<5} # {:5} {:>9n} {:>9n} | {:5} {:>9n} {:>9n} {} | {:>7n} {:.2f} {:.2f} {}\n", si,
			ss[0], atoi(ss[1].c_str()), atoi(ss[2].c_str()), 
			ss[6], atoi(ss[7].c_str()), atoi(ss[8].c_str()), ss[5],
			atoi(ss[17].c_str()), atof(ss[25].c_str()), atof(ss[26].c_str()), 
			ss[16]
		);

		int i_pass = 0, i_total_fails = 0;

		auto aln = align(refa, refb);
		auto alns = aln.trim().max_sum();
		for (auto &a: alns) {
			auto err = a.calculate_error();
			out += fmt::format("      > {:5} {:>9n} {:>9n} | {:5} {:>9n} {:>9n} {} | {:>7n} {:.2f} {:.2f}\n", 
				a.chr_a, a.start_a, a.end_a,
				a.chr_b, a.start_b, a.end_b, ss[5],
				a.alignment.size(),
				err.mis_error() / 100, err.gap_error() / 100
			);

			// we need to padd it properly to assure their equality!
			// eprn("{} {}", a.a.size(), a.b.size());
			// eprn("{} {}", refa.size(), refb.size());
			if (a.a.size() < a.b.size()) {
				a.a += ref[ca].substr(sa + a.end_a, a.b.size() - a.a.size());
			}
			if (a.a.size() > a.b.size()) {
				if (!rb) {
					a.b += ref[cb].substr(sb + a.end_b, a.a.size() - a.b.size());
				} else {
					int real_start = sb + refb.size() - a.end_b;
					string rc = ref[cb].substr(real_start - (a.a.size() - a.b.size()), a.a.size() - a.b.size());
					reverse(rc.begin(), rc.end());
					transform(rc.begin(), rc.end(), rc.begin(), rev_dna);
					a.b += rc;
				}
			}
			assert(a.a.size() == a.b.size());
			Hash ha(a.a);
			Hash hb(a.b);

			bool success = false;

			TREE_t tree;
			vector<Hit> mappings;

			// 1. Try normal Sedef mapping
			for (int i = OFF; i < ha.seq.size(); i += 250) {
				auto m = search(i, ha, hb, tree, true, MIN_READ_SIZE);
				mappings.insert(mappings.end(), m.begin(), m.end());
				for (auto &pp: m) { 
					if (pp.reason.substr(0, 2) == "OK" 
						&& overlap(pp.p, pp.q, OFF, OFF + a.a.size()) >= MIN_ID
						&& overlap(pp.i, pp.j, OFF, OFF + a.b.size()) >= MIN_ID)
					{
						success = true;
						break;
					}
				}
				if (success) break;
			}
			if (success) {
				out += fmt::format("\tExtn:  OK\n");
				i_pass++;
				continue;
			}
			out += fmt::format("\tExtn:  FAIL\n");
			out += print_mappings(mappings, ca, cb, a.a.size(), a.b.size());

			// 2. Try full mappings without extension
			mappings = search(OFF, ha, hb, tree, true, max(a.a.size(), a.b.size()));
			for (auto &pp: mappings) { 
				if (pp.reason.substr(0, 2) == "OK" 
					&& overlap(pp.p, pp.q, OFF, OFF + a.a.size()) >= MIN_ID
					&& overlap(pp.i, pp.j, OFF, OFF + a.b.size()) >= MIN_ID)
				{
					success = true;
					break;
				}
				if (success) break;
			}
			if (success) {
				out += fmt::format("\tFull: OK\n");
				continue;
			}
			out += fmt::format("\tFull: FAIL\n");
			out += a.print();
			out += print_mappings(mappings, ca, cb, a.a.size(), a.b.size());		
			
			i_total_fails++;	
		}

		#pragma omp critical
		{
			eprnn("{}", out);
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
