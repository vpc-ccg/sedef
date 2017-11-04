/// 786

#pragma once

#include <string>
#include "common.h"
#include "fasta.h"

string getfasta(FastaReference &fr, const std::string &chrom, int start, int end, bool rc);
void align(std::string ref_path, std::string bed_path, int resume_after);
