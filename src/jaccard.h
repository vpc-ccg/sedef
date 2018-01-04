/// 786

/******************************************************************************/

#pragma once

/******************************************************************************/

#include <string>

#include "common.h"
#include "hash.h"
#include "search.h"

/******************************************************************************/

void search_main(std::string ref_path, std::string query_path, bool is_complement);

std::string print_mapping(Hit &pp, 
	bool is_complement, 
	const Index &query, 
	const Index &ref);
