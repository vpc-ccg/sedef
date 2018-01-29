/// 786
/// Adapted from SCALCE source code 
/// https://github.com/sfu-compbio/scalce

/******************************************************************************/

#pragma once

/******************************************************************************/

#include <list>
#include <vector>
#include <string>
#include <cmath>

#include "common.h"
#include "hit.h"

using namespace std;

/******************************************************************************/

std::vector<Hit> fast_align(const std::string &query, const std::string &ref, 
	const Hit &orig, int kmer_size = 12);
