/// 786

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

void refine_chains(std::vector<Hit> &anchors, 
	const std::string &qseq, 
	const std::string &rseq,
	const Hit &orig);
