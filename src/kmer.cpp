#include "kmer.h"
#include "util.h"
#include <sstream>
#include <string.h>

Kmer::Kmer(string filename, Options* opt)
{
    mFastaReader = NULL;
    mOptions = opt;
    init(filename);
}

Kmer::~Kmer()
{
    if(mFastaReader) {
        delete mFastaReader;
        mFastaReader = NULL;
    }
}

void Kmer::init(string filename)
{
    mFastaReader = new FastaReader(filename);
    mFastaReader->readAll();

    map<string, string> kmers = mFastaReader->contigs();
    map<string, string>::iterator iter;

    bool initialized = false;
    for(iter = kmers.begin(); iter != kmers.end() ; iter++) {
        string seq = iter->second;

        if(!initialized) {
            initialized = true;
            if(mOptions->kmerKeyLen == 0)
                mOptions->kmerKeyLen = seq.length();
        }
        if(seq.length() != mOptions->kmerKeyLen) {
            cerr << "KMER length must be " << mOptions->kmerKeyLen << ", skipped " << seq << endl;
            continue;
        }
        bool valid = true;
        uint64 kmer64 = seq2uint64(seq, 0, seq.length(), valid);
        if(valid) {
            mKmerHits[kmer64] = 0;
            mNames[kmer64] = iter->first;
            mSequences[kmer64] = iter->second;
        } else {
            cerr << iter->first << ": " << seq << " skipped" << endl;
        }
    }

    if(mKmerHits.size() == 0) {
        error_exit("No unique KMER specified!");
    }
}

void Kmer::report() {
    unordered_map<uint64, uint32>::iterator iter;
    for(iter = mKmerHits.begin(); iter != mKmerHits.end(); iter++) {
        uint64 kmer64 = iter->first;
        cerr << mNames[kmer64] << ": " << mSequences[kmer64] << ": " << iter->second << endl;
    }

    double meanHit = getMeanHit();
    cerr << endl;
    cerr << "Mean coverage: " << meanHit << endl<<endl;
    if(meanHit >= mOptions->positiveThreshold)
        cerr << "Result: POSITIVE";
    else
        cerr << "Result: NEGATIVE";
    cerr << " (" << "threshold: " <<  mOptions->positiveThreshold << ")" << endl;
}

double Kmer::getMeanHit() {
    if(mKmerHits.size() == 0)
        return 0.0;

    double total = 0;
    unordered_map<uint64, uint32>::iterator iter;
    for(iter = mKmerHits.begin(); iter != mKmerHits.end(); iter++) {
        total += iter->second;
    }
    return total / (double) mKmerHits.size();
}

bool Kmer::add(uint64 kmer64) {
    unordered_map<uint64, uint32>::iterator iter = mKmerHits.find(kmer64);
    if(iter != mKmerHits.end()) {
        iter->second++;
        return true;
    }
    return false;
}


string Kmer::getPlotX() {
    stringstream ss;
    unordered_map<uint64, uint32>::iterator iter;
    int first = true;
    for(iter = mKmerHits.begin(); iter != mKmerHits.end(); iter++) {
        if(first) {
            first = false;
        } else 
            ss << ",";

        uint64 kmer64 = iter->first;
        ss << "\"" << mNames[kmer64] << ": " << mSequences[kmer64]<< "\"";
    }
    return ss.str();
}

string Kmer::getPlotY() {
    stringstream ss;
    unordered_map<uint64, uint32>::iterator iter;
    int first = true;
    for(iter = mKmerHits.begin(); iter != mKmerHits.end(); iter++) {
        if(first) {
            first = false;
        } else 
            ss << ",";

        ss << iter->second;
    }
    return ss.str();
}

int Kmer::getKmerCount() {
    return mKmerHits.size();
}

void Kmer::reportJSON(ofstream& ofs) {
    unordered_map<uint64, uint32>::iterator iter;
    int first = true;
    for(iter = mKmerHits.begin(); iter != mKmerHits.end(); iter++) {
        if(first) {
            first = false;
        } else 
            ofs << "," << endl;

        uint64 kmer64 = iter->first;
        ofs << "\t\t\t\"" << mNames[kmer64] << "_" << mSequences[kmer64]<< "\"";
        ofs << ":" << iter->second;
    }
    ofs << endl;
}

uint64 Kmer::seq2uint64(string& seq, uint32 pos, uint32 len, bool& valid) {
    uint64 key = 0;
    for(uint32 i=0; i<len; i++) {
        key = (key << 2);
        switch(seq[pos +i]) {
            case 'A':
                key += 0;
                break;
            case 'T':
                key += 1;
                break;
            case 'C':
                key += 2;
                break;
            case 'G':
                key += 3;
                break;
            case 'N':
            default:
                valid = false;
                return 0;
        }
    }
    valid = true;
    return key;
}