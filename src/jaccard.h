/// 786

/******************************************************************************/

#pragma once

/******************************************************************************/

#include <string>

#include "common.h"
#include "search.h"

/******************************************************************************/

void jaccard_search (std::string ref_path, std::string query_path, bool is_complement);

std::string print_mapping(Hit &pp, 
	bool is_complement, 
	const std::string &query_chr, 
	const std::string &ref_chr,
	int ref_size);