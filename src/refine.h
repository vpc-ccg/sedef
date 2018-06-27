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

void refine_chains(std::vector<Hit> &anchors, 
	const std::string &qseq, 
	const std::string &rseq,
	const Hit &orig);
