#ifndef GENOMES_H
#define GENOMES_H

// includes
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include "common.h"
#include "fastareader.h"
#include <vector>
#include <list>
#include <set>
#include <unordered_map>
#include "options.h"

using namespace std;

class MapResult{

public:
    MapResult(){
        mapped = false;
        start = 0;
        len = 0;
        ed = 0x7FFFFFFF; // initialized with a very large ED
    }
public:
    bool mapped;
    uint32 start;
    uint32 len;
    uint32 ed; // edit distance
};

class Genomes
{
public:
    Genomes(string fastaFile, Options* opt);
    ~Genomes();

    void cover(int id, uint32 pos, uint32 len, uint32 ed, float frac);
    bool hasKey(uint64 key);
    bool align(string& seq);
    void report();
    void reportJSON(ofstream& ofs);
    void reportHtml(ofstream& ofs);

    static uint32 packIdPos(uint32 id, uint32 position);
    static void unpackIdPos(uint32 data,uint32& id, uint32& pos);

private:
    void init();
    void buildKmerTable();
    void addKmer(uint64 key, uint32 id, uint32 pos);
    void initLowComplexityKeys();
    MapResult mapToGenome(string& seq, uint32 seqPos, string& genome, uint32 genomePos);
    void initBloomFilter();
    string getPlotX(int id);
    string getCoverageY(int id);
    string getEditDistanceY(int id);
    void initBinSize();
    double getCoverageRate(int id);

private:
    int mGenomeNum;
    FastaReader* mFastaReader;
    vector<string> mSequences;
    vector<string> mNames;
    vector<vector<float>> mCoverage;
    vector<vector<float>> mEditDistance;
    vector<long> mTotalEditDistance;
    vector<long> mReads;
    vector<long> mBases;
    // unit32 = 8 bits genome id + 24 bits positions
    unordered_map<uint64, list<uint32>> mKmerTable; 
    set<uint64> mLowComplexityKeys;
    Options* mOptions;
    long mHitCount;
    long mMissedCount;
    char* mBloomFilterArray;
};


#endif