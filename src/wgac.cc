/// 786

/******************************************************************************/

#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <thread>
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

double overlap(int sa, int ea, int sb, int eb)
{
	int o = min(ea, eb) - max(sa, sb);
	// eprn("overlap: {} {} / {} {} = {}", sa, ea, sb, eb, o);
	o = max(0, o);

	return o / double(eb - sb);
}

extern unique_ptr<AHOAutomata> aho;

void check_wgac(string bed_path) 
{
   aho = make_unique<AHOAutomata>();

	string s;

	eprnn("Loading reference... ");
	unordered_map<string, string> ref;
	ifstream ff("data/hg19/chr1.fa");
	string chr;
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
	vector<string> lines;
	while (getline(fin, s)) {
		lines.push_back(s);
	}

	int total = 0, pass = 0, total_fails = 0;
	auto xx =  set<int>{152, 1036 };
	#pragma omp parallel for
	for (int si = 0; si < lines.size(); si++) {
		// if (xx.find(si) == xx.end()) continue;

		auto ss = split(lines[si], '\t');

		string ca = ss[0];
		int sa = atoi(ss[1].c_str()), ea = atoi(ss[2].c_str());
		string cb = ss[6];
		int sb = atoi(ss[7].c_str()), eb = atoi(ss[8].c_str());
		bool rb = (ss[5][0] == '_');

		if (ca != "chr1" || cb != "chr1") continue;

		const int OFF = 0;
		int max_to_extract = max(ea - sa, eb - sb);
		auto refa = ref[ca].substr(sa - OFF, max_to_extract + 2 * OFF);
		auto refb = ref[cb].substr(sb - OFF, max_to_extract + 2 * OFF);
		if (rb) {
			reverse(refb.begin(), refb.end());
        	transform(refb.begin(), refb.end(), refb.begin(), rev_dna);
		}

		Hash ha(refa);
		Hash hb(refb);

		int passes = 0;

		string out;
		out += fmt::format("{:<5} # {:5} {:>9} {:>9} | {:5} {:>9} {:>9} {} | {:>7} {} {} {}\n", si,
			ss[0], ss[1], ss[2], 
			ss[6], ss[7], ss[8], ss[5],
			ss[17], ss[16], ss[25], ss[26]);
		
		bool success = false;
		TREE_t tree;
		auto mapping = search(OFF, ha, hb, tree, true, max_to_extract);
		for (auto &pp: mapping) { 
			if (pp.reason.substr(0, 2) == "OK" 
				&& overlap(pp.p, pp.q, OFF, OFF + ea - sa) >= .9
				&& overlap(pp.i, pp.j, OFF, OFF + eb - sb) >= .9)
			{
				success = true;
				break;
			}
			if (success) break;
		}
		if (!success) {
			out += fmt::format("\tFull: FAIL\n");

			sort(mapping.begin(), mapping.end(), [&](Hit pp, Hit pq) {
				auto ta = make_pair(overlap(pp.p, pp.q, OFF, OFF + ea - sa), overlap(pp.i, pp.j, OFF, OFF + eb - sb));
				auto tb = make_pair(overlap(pq.p, pq.q, OFF, OFF + ea - sa), overlap(pq.i, pq.j, OFF, OFF + eb - sb));
				return ta > tb;
			});
			int cnt = count_if(mapping.begin(), mapping.end(), [](Hit pp){ return pp.reason.substr(0, 2) == "OK"; });
			int prn = 0;
			
			double prev_score = -1;
			for (auto &pp: mapping) { 
				if (!cnt || pp.reason.substr(0, 2) == "OK") {
					double score = overlap(pp.p, pp.q, OFF, OFF + ea - sa) + overlap(pp.i, pp.j, OFF, OFF + eb - sb);
					if (prev_score != -1 && prev_score / score > 2) break;
					out += "\t\t" + fmt::format("{:.2f} {:.2f} ~ ", 
						overlap(pp.p, pp.q, OFF, OFF + ea - sa),
						overlap(pp.i, pp.j, OFF, OFF + eb - sb)) + print_mapping(pp, 0, ca, cb, hb); 
					prev_score = score;
				}
			}

			// eprnn("{}", out);
			auto aln = align(refa, refb);
			out += aln.print();
			auto alns = aln.trim().max_sum();
			int tt(0);
			for (auto &a: alns) {
				out += fmt::format("   >> Trim {}:\n", ++tt);
				out += a.print();
			}
			// exit(0);
		} else {
			passes += 10;
			out += fmt::format("\tFull: OK\n");
		}
		
		success = false;
		vector<Hit> mappings;
		tree.clear();
		for (int i = OFF; i < refa.size(); i += 250) {
			mapping = search(i, ha, hb, tree, true, 1000);
			mappings.insert(mappings.end(), mapping.begin(), mapping.end());
			for (auto &pp: mapping) { 
				if (pp.reason.substr(0, 2) == "OK" 
					&& overlap(pp.p, pp.q, OFF, OFF + ea - sa) >= .9
					&& overlap(pp.i, pp.j, OFF, OFF + eb - sb) >= .9)
				{
					success = true;
					break;
				}
			}
			if (success) break;
		}
		if (!success) {
			mapping = mappings;
			out += fmt::format("\tExt:  FAIL\n");

			sort(mapping.begin(), mapping.end(), [&](Hit pp, Hit pq) {
				auto ta = make_pair(overlap(pp.p, pp.q, OFF, OFF + ea - sa), overlap(pp.i, pp.j, OFF, OFF + eb - sb));
				auto tb = make_pair(overlap(pq.p, pq.q, OFF, OFF + ea - sa), overlap(pq.i, pq.j, OFF, OFF + eb - sb));
				return ta > tb;
			});
			int cnt = count_if(mapping.begin(), mapping.end(), [](Hit pp){ return pp.reason.substr(0, 2) == "OK"; });
			int prn = 0;
			
			double prev_score = -1;
			for (auto &pp: mapping) { 
				if (!cnt || pp.reason.substr(0, 2) == "OK") {
					double score = overlap(pp.p, pp.q, OFF, OFF + ea - sa) + overlap(pp.i, pp.j, OFF, OFF + eb - sb);
					if (prev_score != -1 && prev_score / score > 2) break;
					out += "\t\t" + fmt::format("{:.2f} {:.2f} ~ ", 
						overlap(pp.p, pp.q, OFF, OFF + ea - sa),
						overlap(pp.i, pp.j, OFF, OFF + eb - sb)) + print_mapping(pp, 0, ca, cb, hb); 
					prev_score = score;
				}
			}
		} else {
			passes++;
			out += fmt::format("\tExt:  OK\n");
		}

		#pragma omp critical
		{
			if (passes != 11) {
				prnn("{}", out);
				// break;
			} else {
				pass++;
			}
			if (passes < 10) total_fails++;
			total++;
		}
	}

	eprn("passed {}/{} total_fail {}", pass, total, total_fails);
}
