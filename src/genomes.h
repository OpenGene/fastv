#ifndef GENOMES_H
#define GENOMES_H

// includes
#include "common.h"
#include "fastareader.h"
#include <vector>
#include <list>
#include "options.h"

using namespace std;

class Genomes
{
public:
    Genomes(string fastaFile, Options* opt);
    ~Genomes();

    void cover(int id, uint32 pos, uint32 len);

    static uint32 packIdPos(uint32 id, uint32 position);
    static void unpackIdPos(uint32 data,uint32& id, uint32& pos);

private:
    void init();
    void buildKmerTable();
    void addKmer(uint64 key, uint32 id, uint32 pos);

private:
    int mGenomeNum;
    FastaReader* mFastaReader;
    vector<string> mSequences;
    vector<string> mNames;
    vector<vector<uint16>> mCoverage;
    // unit32 = 10 bits genome id + 22 bits positions
    map<uint64, list<uint32>> mKmerTable; 
    Options* mOptions;
};


#endif