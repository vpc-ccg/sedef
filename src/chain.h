/// 786

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Author: inumanag

/******************************************************************************/

#pragma once

/******************************************************************************/

#include <list>
#include <vector>
#include <string>
#include <cmath>

#include "common.h"
#include "hit.h"

/******************************************************************************/

std::vector<Hit> fast_align(const std::string &query, const std::string &ref, 
	const Hit &orig, int kmer_size);
