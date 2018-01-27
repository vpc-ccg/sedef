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

vector<Hit> merge(vector<Hit> &hits, const int merge_dist = 100) {
	vector<Hit> results;

	sort(hits.begin(), hits.end(), [](const Hit &a, const Hit &b) {
		return 
			tie(a.ref->is_rc, a.query->name, a.query_start, a.ref->name, a.ref_start) <
			tie(b.ref->is_rc, b.query->name, b.query_start, b.ref->name, b.ref_start);
	});

	Hit rec, prev;
	int wcount = 0;
	size_t len = 0;
	ssize_t nread;
	multimap<int, Hit> windows;
	for (auto &rec: hits) {
		assert(!rec.query->is_rc);
		if (rec.query->name == rec.ref->name && rec.query_start == rec.ref_start && rec.query_end == rec.ref_end &&
			rec.query->is_rc == rec.ref->is_rc)
			continue;
		if ((&rec - &hits[0]) == 0) {
			windows.emplace(rec.ref_end, rec);
			prev = rec;
			wcount++;
		} else if (prev.query_end + merge_dist < rec.query_start || prev.query->name != rec.query->name || 
				prev.ref->name != rec.ref->name || prev.ref->is_rc != rec.ref->is_rc) 
		{
			for (auto it: windows)
				results.push_back(it.second);
			windows.clear();
			windows.emplace(rec.ref_end, rec);
			prev = rec;
			wcount++;
		} else {
			bool needUpdate = 1;
			while (needUpdate) {
				auto start_loc = windows.lower_bound(rec.ref_start - merge_dist);
				needUpdate = 0;
				while (start_loc != windows.end()) {
					if (start_loc->second.query_end + merge_dist < rec.query_start ||
							  start_loc->second.ref_end < rec.ref_start - merge_dist ||
							  start_loc->second.ref_start > rec.ref_end + merge_dist) 
					{
						 start_loc++;
						 continue;
					}
					needUpdate = 1;
					rec.query_end = max(rec.query_end, start_loc->second.query_end);
					rec.ref_end = max(rec.ref_end, start_loc->second.ref_end);
					rec.query_start = min(rec.query_start, start_loc->second.query_start);
					rec.ref_start = min(rec.ref_start, start_loc->second.ref_start);
					windows.erase(start_loc++);
				}
			}
			windows.emplace(rec.ref_end, rec);
		}
		rec.query_end = max(rec.query_end, prev.query_end);
		prev = rec;
	}
	for (auto it: windows)
		results.push_back(it.second);
	return results;
}

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
				int w = max(h.query_end - h.query_start, h.ref_end - h.ref_start);
				w = min(15000, 5 * w);
				h.query_start = max(0, h.query_start - w);
				h.query_end += w;
				h.ref_start = max(0, h.ref_start - w);
				h.ref_end += w;
			}

			hits.push_back(h);
			nhits++;
		}
		eprn("Read {} alignments in {}", nhits, file);
	}

	eprn("Read total {} alignments", hits.size());
	if (extend) {
		hits = merge(hits);
		eprn("After merging remaining {} alignments", hits.size());
	}

	vector<vector<Hit>> bins(10000);
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

