/// 786

#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include "fmt/fmt/format.h"
#include "fmt/fmt/format.cc"
#include "ksw2.h"
#include "ksw2_extz2_sse.c"
using namespace std;
#define prn(f, ...)    fmt::print(f "\n", __VA_ARGS__)

static char dna_lookup[128] = {
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
	4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
	3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 
	4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

static char rev_dna_lookup[128] = {
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
	4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 2, 4, 4, 4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
	0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 2, 4, 4, 4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 
	4, 4, 4, 4, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
}

struct FastaIndexEntry {
	FastaIndexEntry(string name, int length, long long offset, int line_blen, int line_len);
    string name;  // sequence name
    int length;  // length of sequence
    long long offset;  // bytes offset of sequence from start of file
    int line_blen;  // line length in bytes, sequence characters
    int line_len;  // line length including newline
};

FastaIndexEntry::FastaIndexEntry(string name, int length, long long offset, int line_blen, int line_len)
    : name(name)
    , length(length)
    , offset(offset)
    , line_blen(line_blen)
    , line_len(line_len)
{}


struct FastaIndex : public map<string, FastaIndexEntry> {
    FastaIndex() {}
    ~FastaIndex(void);
    vector<string> sequenceNames;
    void readIndexFile(const string &fname);
    ifstream indexFile;
    FastaIndexEntry entry(string key);
};

FastaIndex::~FastaIndex(void) {
    indexFile.close();
}

void FastaIndex::readIndexFile(const string &fname) {
    string line;
    long long linenum = 0;
    indexFile.open(fname.c_str(), ifstream::in);
    if (indexFile.is_open()) {
        while (getline (indexFile, line)) {
            ++linenum;
            // the fai format defined in samtools is tab-delimited, every line being:
            // fai->name[i], (int)x.len, (long long)x.offset, (int)x.line_blen, (int)x.line_len
            vector<string> fields = split(line, '\t');
            if (fields.size() == 5) {  // if we don't get enough fields then there is a problem with the file
                // note that fields[0] is the sequence name
                char* end;
                string name = split(fields[0], ' ').at(0);  // key by first token of name
                sequenceNames.push_back(name);
                this->insert(make_pair(name, FastaIndexEntry(
                				fields[0], atoi(fields[1].c_str()),
							    strtoll(fields[2].c_str(), &end, 10),
							    atoi(fields[3].c_str()),
							    atoi(fields[4].c_str())) ));
            } else {
                cerr << "Warning: malformed fasta index file " << fname <<  "does not have enough fields @ line " << linenum << endl;
                exit(1);
            }
        }
    } else {
        cerr << "could not open index file " << fname << endl;
        exit(1);
    }
}

FastaIndexEntry FastaIndex::entry(string name) {
    FastaIndex::iterator e = this->find(name);
    if (e == this->end()) {
        cerr << "unable to find FASTA index entry for '" << name << "'" << endl;
        exit(1);
    } else {
        return e->second;
    }
}

struct FastaReference {
	FastaReference(string reffilename);
	~FastaReference();
    void open(string reffilename, bool usemmap = false,
		  bool useFullHeader = false);
    string filename;
    FILE* file;
    void* filemm;
    size_t filesize;
    FastaIndex* index;
    // potentially useful for performance, investigate
    string getSubSequence(string seqname, int start, int length);
};

FastaReference::FastaReference(string reffilename) {
    filename = reffilename;
    if (!(file = fopen(filename.c_str(), "r"))) {
        cerr << "could not open " << filename << endl;
        exit(1);
    }
    index = new FastaIndex();
    struct stat stFileInfo; 
    string indexFileName = filename + ".fai"; 
    // if we can find an index file, use it
    if(stat(indexFileName.c_str(), &stFileInfo) == 0) {
        // check if the index file is older than the FASTA file
        struct stat index_attrib, fasta_attrib;
        stat(indexFileName.c_str(), &index_attrib);
        stat(reffilename.c_str(), &fasta_attrib);
        if (fasta_attrib.st_mtime > index_attrib.st_mtime) {
            cerr << "Warning: the index file is older than the FASTA file." << endl;
        }
        index->readIndexFile(indexFileName);
    } 
    int fd = fileno(file);
    struct stat sb;
    if (fstat(fd, &sb) == -1)
        cerr << "could not stat file" << filename << endl;
    filesize = sb.st_size;
    // map the whole file
    filemm = mmap(NULL, filesize, PROT_READ, MAP_SHARED, fd, 0);
}

FastaReference::~FastaReference(void) {
    fclose(file);
    munmap(filemm, filesize);
    delete index;
}

string FastaReference::getSubSequence(string seqname, int start, int length) {
    FastaIndexEntry entry = index->entry(seqname);
    if (start < 0 || length < 1) {
        cerr << "Error: cannot construct subsequence with negative offset or length < 1" << endl;
        exit(1);
    }
    // we have to handle newlines
    // approach: count newlines before start
    //           count newlines by end of read
    //             subtracting newlines before start find count of embedded newlines
    if (start + length > entry.length) {
    	length = entry.length - start;
    }
    int newlines_before = start > 0 ? (start - 1) / entry.line_blen : 0;
    int newlines_by_end = (start + length - 1) / entry.line_blen;
    int newlines_inside = newlines_by_end - newlines_before;
    int seqlen = length + newlines_inside;
    char* seq = (char*) calloc (seqlen + 1, sizeof(char));

    memcpy(seq, (char*) filemm + entry.offset + newlines_before + start, seqlen);
    seq[seqlen] = '\0';
    char* pbegin = seq;
    char* pend = seq + (seqlen/sizeof(char));
    pend = remove(pbegin, pend, '\n');
    pend = remove(pbegin, pend, '\0');
    string s = seq;
    free(seq);
    s.resize((pend - pbegin)/sizeof(char));
    return s;
}

string getfasta(FastaReference &fr, const string &chrom, int start, int end, bool rc) 
{
	/* Make sure that we can open all of the files successfully*/
	string s = fr.getSubSequence(chrom, start, end - start);
	// prn("{} // {}", s, s.size());
	char *plook = dna_lookup;
	if (rc) {
		reverse(s.begin(), s.end());
		plook = rev_dna_lookup;
	}
	// prn("{} // {}", s, s.size());
	// string q;
    for (int i = 0; i < s.size(); i++) {
    	s[i] = plook[s[i]];
    	// q+="ACGTN"[s[i]];
    }
	// prn("{}\n{}", rc, q);
    return s;
}

string align(string &tseq, string &qseq, int sc_mch, int sc_mis, int gapo, int gape)
{
	int8_t a = (int8_t)sc_mch, b = sc_mis < 0 ? (int8_t)sc_mis : (int8_t)(-sc_mis); // a>0 and b<0
	int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
	string cigar;
	const int STEP = 50000;
	for (int SP = 0; SP < min(tseq.size(), qseq.size()); SP += STEP) {
		ksw_extz_t ez;
		ksw_extz2_sse(0, 
			min(STEP, (int)(qseq.size() - SP)), (const uint8_t*)(qseq.c_str() + SP), 
			min(STEP, (int)(tseq.size() - SP)), (const uint8_t*)(tseq.c_str() + SP), 
			5, mat, // M; MxM matrix
			gapo, gape, 
			-1, -1, // band width; off-diagonal drop-off to stop extension (-1 to disable)
			0, &ez);
		
		for (int i = 0; i < ez.n_cigar; ++i) 
			cigar += fmt::format("{}{}", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
		free(ez.cigar);
		cigar += ";";
	}
	return cigar;
}


std::chrono::time_point<std::chrono::high_resolution_clock> t_start;
inline int TIME()
{
    auto end = chrono::high_resolution_clock::now();
    return chrono::duration_cast<chrono::seconds>(end - t_start).count();
}


void align(FastaReference &fr, int line, const string &s)
{
	vector<string> ss = split(s, '\t');

	char ca[500]; int sa, ea;
	char cb[500]; int sb, eb;
	char sta[10], stb[10];
	sscanf(s.c_str(), "%s %d %d %s %d %d %*s %s %s", ca, &sa, &ea, cb, &sb, &eb, sta, stb);
	assert(sta[0] == '+');
	bool rc = (stb[0] == '-');


	string fa = getfasta(fr, string(ca), sa - 1, ea - 1, false);
	string fb = getfasta(fr, string(cb), sb - 1, eb - 1, rc);

	string aln = align(fa, fb, 5, -4, 40, 1);

	cerr << "\t" << TIME() << "s\t" << line << endl;
	cout << s << "\t" << aln << endl;
}

//389,419,110/20000000
// 19.47

int main(int argc, char **argv)  // megaalign fasta/hg19.fa output/search/buckets 500 20000000 0
{ // usage: reference.fa bucket.bed
	std::ios_base::sync_with_stdio(0);

	system(">&2 date");

	if (argc < 3) exit(1);
	t_start = chrono::high_resolution_clock::now();

	FastaReference fr(argv[1]);
	ifstream fin(argv[2]);

    int startfrom = -1;
    if (argc > 3) startfrom = atoi(argv[3]);

	int nl = 0;
	string s;
	while (getline(fin, s)) {
		if (nl > startfrom) align(fr, nl, s);
        nl++;
	}

	cerr << fmt::format("Finished bcuket {}: read {} lines", string(argv[2]), nl) << endl;
	cerr << "Done!" << endl;

	system(">&2 date");

	return 0;
}
