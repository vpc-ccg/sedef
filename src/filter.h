/// 786

/******************************************************************************/

#pragma once

/******************************************************************************/

#include <list>
#include <vector>
#include <string>
#include <cmath>

#include "common.h"
#include "search.h"

using namespace std;

/******************************************************************************/

pair<bool, string> filter(
	const string &q, int q_pos, int q_len, 
	const string &r, int r_pos, int r_len);

