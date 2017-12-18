/// 786
/// Taken from Bedtools source code 

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
	string getSubSequence(string seqname, int start, int length);
};
