#ifndef KMER_H
#define KMER_H

// includes
#include "common.h"
#include <vector>
#include <map>
#include "fastareader.h"
#include "options.h"

using namespace std;

class Kmer
{
public:
    Kmer(string filename, Options* opt);
    ~Kmer();
    void init(string filename);
    bool add(uint64 kmer64);
    void report();
    double getMeanHit();

    static uint64 seq2uint64(string& seq, uint32 pos, uint32 len, bool& valid);

private:
    map<uint64, uint32> mKmerHits;
    FastaReader* mFastaReader;
    map<uint64, string> mNames;
    map<uint64, string> mSequences;
    Options* mOptions;
};


#endif