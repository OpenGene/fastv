#include "kmer.h"
#include "util.h"

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

    for(iter = kmers.begin(); iter != kmers.end() ; iter++) {
        string seq = iter->second;
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
}

void Kmer::report() {
    map<uint64, uint32>::iterator iter;
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
    map<uint64, uint32>::iterator iter;
    for(iter = mKmerHits.begin(); iter != mKmerHits.end(); iter++) {
        total += iter->second;
    }
    return total / (double) mKmerHits.size();
}

bool Kmer::add(uint64 kmer64) {
    map<uint64, uint32>::iterator iter = mKmerHits.find(kmer64);
    if(iter != mKmerHits.end()) {
        iter->second++;
        return true;
    }
    return false;
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