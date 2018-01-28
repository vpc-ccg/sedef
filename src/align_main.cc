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
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "align.h"
#include "align_main.h"
#include "common.h"
#include "fasta.h"
#include "chain.h"
#include "hit.h"
#include "merge.h"

using namespace std;

/******************************************************************************/

const double EXTEND_RATIO = 5;
const int MAX_EXTEND = 15000;
const int MERGE_DIST = 250;

/******************************************************************************/

auto stat_file(const string &path)
{
	struct stat path_stat;
	int s = stat(path.c_str(), &path_stat);
	assert(s == 0);
	return path_stat.st_mode;
}

auto bucket_alignments(const string &bed_path, int nbins, string output_dir = "", bool extend=false)
{
	vector<string> files;
	if (S_ISREG(stat_file(bed_path))) {
		files.push_back(bed_path);
	} else if (S_ISDIR(stat_file(bed_path))) {
		glob_t glob_result;
		glob((bed_path + "/*.bed").c_str(), GLOB_TILDE, NULL, &glob_result);
		for (int i = 0; i < glob_result.gl_pathc; i++) {
			string f = glob_result.gl_pathv[i];
			if (S_ISREG(stat_file(f))) {
				files.push_back(f);
			}
		}
	} else {
		throw fmt::format("Path {} is neither file nor directory", bed_path);
	}

	vector<Hit> hits;
	for (auto &file: files) {
		ifstream fin(file.c_str());
		if (!fin.is_open()) {
			throw fmt::format("BED file {} does not exist", bed_path);
		}

		int nhits = 0;
		string s;
		while (getline(fin, s)) {
			Hit h = Hit::from_bed(s);
			if (extend) {
				h.extend(EXTEND_RATIO, MAX_EXTEND);
			}
			hits.push_back(h);
			nhits++;
		}
		eprn("Read {} alignments in {}", nhits, file);
	}

	eprn("Read total {} alignments", hits.size());
	if (extend) {
		hits = merge(hits, MERGE_DIST);
		eprn("After merging remaining {} alignments", hits.size());
	}
	int max_complexity = 0;
	for (auto &h: hits) {
		max_complexity = max(
			max_complexity, 
			(int)sqrt(double(h.query_end - h.query_start) * double(h.ref_end - h.ref_start))
		);
	}
	vector<vector<Hit>> bins(max_complexity / 1000 + 1);
	for (auto &h: hits) {
		int complexity = sqrt(double(h.query_end - h.query_start) * double(h.ref_end - h.ref_start));
		assert(complexity / 1000 < bins.size());
		bins[complexity / 1000].push_back(h);
	}

	vector<vector<Hit>> results(nbins);
	int bc = 0;	
	for (auto &bin: bins) {
		for (auto &hit: bin) {
			results[bc].push_back(hit);
			bc = (bc + 1) % nbins;
		}
	}

	if (output_dir != "") {
		int count = 0;
		for (auto &bin: results) {
			string of = output_dir + fmt::format("/bucket_{:04d}", count++);
			ofstream fout(of.c_str());
			if (!fout.is_open()) {
				throw fmt::format("Cannot open file {} for writing", of);
			}
			for (auto &h: bin) {
				fout << h.to_bed(false) << endl;
			}
			fout.close();
			eprn("Wrote {} alignments in {}", bin.size(), of);
		}
	}

	return results;
}

void generate_alignments(const string &ref_path, const string &bed_path, int kmer_size) 
{
	auto T = cur_time();

	auto schedule = bucket_alignments(bed_path, 1);
	FastaReference fr(ref_path);

	int lines = 0, total = 0;
	for (auto &s: schedule) 
		total += s.size();

	eprn("Using k-mer size {}", kmer_size);

	int total_written = 0;
	for (int i = 0; i < schedule.size(); i++) {
		for (auto &h: schedule[i]) {
			string fa = fr.get_sequence(h.query->name, h.query_start, &h.query_end);
			string fb = fr.get_sequence(h.ref->name, h.ref_start, &h.ref_end);
			if (h.ref->is_rc) 
				fb = rc(fb);

			auto alns = fast_align(fa, fb, kmer_size);
			lines++;
			for (auto &hh: alns) {
				hh.query_start += h.query_start;
				hh.query_end += h.query_start;
				if (h.ref->is_rc) {
					swap(hh.ref_start, hh.ref_end);
					hh.ref_start = h.ref_end - hh.ref_start;
					hh.ref_end = h.ref_end - hh.ref_end;
					hh.ref->is_rc = true;
				} else {
					hh.ref_start += h.ref_start;
					hh.ref_end += h.ref_start;
				}
				hh.query->name = h.query->name;
				hh.ref->name = h.ref->name;
				hh.ref->name = h.ref->name;
				total_written++;
				prn("{}", hh.to_bed(false));
			}
			eprnn("\r {} out of {} ({:.1f}, len {}..{})      ", lines, total, pct(lines, total),
				fa.size(), fb.size());
		}
	}

	eprn("\nFinished BED {} in {}s ({} lines, generated {} hits)", bed_path, elapsed(T), lines, total_written);
}

/******************************************************************************/

void align_main(int argc, char **argv)
{
	if (argc < 3) {
		throw fmt::format("Not enough arguments to align");
	}

	string command = argv[0];
	if (command == "bucket") {
		if (argc < 4) {
			throw fmt::format("Not enough arguments to align-bucket");
		}
		bucket_alignments(argv[1], atoi(argv[3]), argv[2], true);
	} else if (command == "generate") {
		generate_alignments(argv[1], argv[2], atoi(argv[3]));
	// } else if (command == "process") {
	// 	postprocess(argv[1], argv[2]);
	} else {
		throw fmt::format("Unknown align command");
	}
}

