#ifndef ALLKMER_H
#define ALLKMER_H

// includes
#include "common.h"
#include <vector>
#include <unordered_map>
#include <map>
#include "fastareader.h"
#include "options.h"
#include "zlib/zlib.h"
#include "common.h"
#include <iostream>
#include <fstream>
#include <mutex>

#define  MTX_COUNT 100
#define COLLISION_FLAG 0xFFFFFFFF

using namespace std;

class KCResult {
public:
    string mName;
    uint64  mHit;
    int mMedianHit;
    double mMeanHit;
    double mCoverage;
    int mKmerCount;
};

class KCHit {
public:
    uint64 mKey64;
    uint32 mID;
    uint32 mHit;
};

class KmerCollection
{
public:
    KmerCollection(string filename, Options* opt);
    ~KmerCollection();
    void init();
    void report();
    void reportJSON(ofstream& ofs);
    void reportHTML(ofstream& ofs);
    bool add(uint64 kmer64);

    uint32 packIdCount(uint32 id, uint32 count);
    void unpackIdCount(uint32 data,uint32& id, uint32& count);
    void stat();

private:
    bool getLine(char* line, int maxLine);
    uint64 makeHash(uint64 key);
    bool eof();
    void makeBitAndMask();
    bool isHighConfidence(KCResult kcr);
private:
    Options* mOptions;
    vector<string> mNames;
    vector<uint64> mHits;
    vector<int> mMedianHits;
    vector<double> mMeanHits;
    vector<double> mCoverage;
    vector<int> mKmerCounts;
    vector<KCResult> mResults;
    int mNumber;
    uint32 mUniqueHashNum;
    uint32* mHashKCH;
    KCHit* mKCHits;
    string mFilename;
    gzFile mZipFile;
    ifstream mFile;
    bool mZipped;
    int mIdBits;
    uint32 mIdMask;
    uint32 mCountMax;
    bool mStatDone;
    uint32 mUniqueNumber;
};


#endif