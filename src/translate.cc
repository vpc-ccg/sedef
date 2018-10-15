/// 786 

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Author: inumanag

/******************************************************************************/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>

#include "common.h"
#include "fasta.h"
#include "extern/argh.h"

using namespace std;

/******************************************************************************/

void translate(const string &ref_path, const string &out_path, const size_t max_size)
{
   auto T = cur_time();

   eprn("Translating {} to {}", ref_path, out_path);

   FastaReference fr(ref_path);

   vector<pair<size_t, string>> vv;
   for (const auto &rf: fr.index) {
      vv.push_back({rf.second.length, rf.second.name});
   }
   sort(vv.begin(), vv.end(), std::greater<pair<size_t, string>>());
      
   FILE *fo = fopen(out_path.c_str(), "w");
   FILE *foi = fopen(fmt::format("{}.sedef-translate", out_path).c_str(), "w");
   int chr_id = 1;
   string chr = "";
   for (auto &v: vv) {
      fprintf(foi, "%s\ttrans%d\t%zu\n", v.second.c_str(), chr_id, chr.size());
      string ref = fr.get_sequence(v.second);      
      chr += ref;

      if (chr.size() > max_size) {
         eprn("Chromosome {:3}, size {:12n}", chr_id, chr.size());
         fprintf(fo, ">trans%d\n", chr_id++);
         for (int i = 0; i < chr.size(); i += 100)
            fprintf(fo, "%s\n", chr.substr(i, 100).c_str());
         chr = "";
      } else {
         chr += string(Globals::Stats::MIN_ASSEMBLY_GAP_SIZE + 1, 'N');
      }
   }
   if (chr.size()) {
      eprn("Chromosome {:3}, size {:12n}", chr_id, chr.size());
      fprintf(fo, ">trans%d\n", chr_id);
      for (int i = 0; i < chr.size(); i += 100)
         fprintf(fo, "%s\n", chr.substr(i, 100).c_str());
   }
   fclose(fo);
   fclose(foi);

   eprn("Done with {} chromosomes", chr_id);   
   eprn("Translating took {:.1f}s", elapsed(T)), T = cur_time();
}

/******************************************************************************/

void trans_main(int argc, char **argv)
{
   argh::parser cmdl;
   cmdl.add_params({"-m", "--max-size"});
   cmdl.parse(argc, argv);
   
   if (!cmdl(1)) {
      throw fmt::format("Not enough arguments to translate");
   }
   
   size_t max_size;
   cmdl({"-m", "--max-size"}, 100 * MB) >> max_size;
   
   translate(cmdl[0], cmdl[1], max_size);
}


