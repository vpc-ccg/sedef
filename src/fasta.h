/// 786
/// Adapted from Bedtools source code 
/// https://github.com/arq5x/bedtools/tree/master/src/fastaFromBed

/******************************************************************************/

#pragma once

/******************************************************************************/

#include <string>
#include <vector>
#include <map>
#include <cstdio>

using namespace std;

/******************************************************************************/

struct FastaIndexEntry {
	string name;
	int length;
	long long offset;
	int line_blen, line_len;

	FastaIndexEntry(string name, int length, long long offset, int line_blen, int line_len): 
		name(name), length(length), offset(offset), 
		line_blen(line_blen), line_len(line_len) 
	{
	}
};

/******************************************************************************/

struct FastaIndex: public map<string, FastaIndexEntry> {
	vector<string> sequenceNames;

	FastaIndex() {};
	FastaIndex(const string &fname);
	FastaIndexEntry entry(const string &key);
};

/******************************************************************************/

struct FastaReference {
	FILE* file;
	void* filemm;
	size_t filesize;
	FastaIndex index;
	
	FastaReference(string filename);
	~FastaReference();
	string get_sequence(string seqname, int start = 0, int end = -1);
};
