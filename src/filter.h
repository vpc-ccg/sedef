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
#include "search.h"

/******************************************************************************/

std::pair<bool, std::string> filter(
	const std::string &q, int q_pos, int q_end, 
	const std::string &r, int r_pos, int r_end);

