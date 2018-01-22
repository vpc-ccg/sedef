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

using namespace std;

/******************************************************************************/

auto stat_file(const string &path)
{
	struct stat path_stat;
	int s = stat(path.c_str(), &path_stat);
	assert(s == 0);
	return path_stat.st_mode;
}

auto bucket_alignments(const string &bed_path, int nbins, string output_dir = "")
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

	vector<vector<Hit>> bins(2000);
	for (auto &file: files) {
		ifstream fin(file.c_str());
		if (!fin.is_open()) {
			throw fmt::format("BED file {} does not exist", bed_path);
		}

		int nhits = 0;
		string s;
		while (getline(fin, s)) {
			Hit h = Hit::from_bed(s);
			int complexity = sqrt(double(h.query_end - h.query_start) * double(h.ref_end - h.ref_start));
			assert(complexity / 1000 < bins.size());
			bins[complexity / 1000].push_back(h);
			nhits++;
		}
		eprn("Read {} alignments in {}", nhits, file);
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

void generate_alignments(const string &ref_path, const string &bed_path) 
{
	auto T = cur_time();

	auto schedule = bucket_alignments(bed_path, 4);
	FastaReference fr(ref_path);

	int lines = 0, total = 0;
	for (auto &s: schedule) 
		total += s.size();

	extern int DEBUG;
	DEBUG = 0; int WW=0;
	#pragma omp parallel for
	for (int i = 0; i < schedule.size(); i++) {
		// auto &h = schedule[i].back();
		for (auto &h: schedule[i]) {
			// if (fmt::format("{} {} {} {}", h.query_start, h.query_end, h.ref_start, h.ref_end) != "62038352 62039402 192191615 192193123") 
				// continue;
			string fa, fb;
			#pragma omp critical 
			{
				fa = fr.get_sequence(h.query->name, h.query_start, h.query_end);
				fb = fr.get_sequence(h.ref->name, h.ref_start, h.ref_end);
			}
			if (h.ref->is_rc) 
				fb = rc(fb);

			auto alns = fast_align(fa, fb);
			#pragma omp critical
			{
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
					WW++;
					prn("{}", hh.to_bed(false));
				}
				eprnn("\r {} out of {} ({:.1f}, len {}..{})", lines, total, pct(lines, total),
					fa.size(), fb.size());
			}
			// if (WW>5) break;
		}
	}

	eprn("\nFinished BED {} in {}s ({} lines)", bed_path, elapsed(T), lines);
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
		bucket_alignments(argv[1], atoi(argv[3]), argv[2]);
	} else if (command == "generate") {
		generate_alignments(argv[1], argv[2]);
	// } else if (command == "process") {
	// 	postprocess(argv[1], argv[2]);
	} else {
		throw fmt::format("Unknown align command");
	}
}

