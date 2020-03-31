#include "virusdetector.h"

VirusDetector::VirusDetector(Options* opt){
    mOptions = opt;

    mKmer = new Kmer(mOptions->kmerFile, opt);
    mGenomes = new Genomes(mOptions->genomeFile, opt);
    mHits = 0;
}

VirusDetector::~VirusDetector(){
    if(mKmer) {
        delete mKmer;
        mKmer = NULL;
    }
    if(mGenomes) {
        delete mGenomes;
        mGenomes = NULL;
    }
}

void VirusDetector::report() {
    if(mKmer) {
        mKmer->report();
    }
}

bool VirusDetector::detect(Read* r) {
    string& seq = r->mSeq.mStr;

    int hitCount = 0;

    int keylen = mOptions->kmerKeyLen;
    int blankBits = 64 - 2*keylen;

    bool valid = true;

    uint32 start = 0;
    uint64 key = Kmer::seq2uint64(seq, start, keylen-1, valid);
    while(valid == false) {
        start++;
        key = Kmer::seq2uint64(seq, start, keylen-1, valid);
        // reach the tail
        if(start > seq.length() - keylen)
            return false;
    }
    for(uint32 pos = start; pos < seq.length() - keylen; pos++) {
        key = (key << 2);
        switch(seq[pos + keylen-1]) {
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
                // we have to skip the segments covering this N
                if(pos >= seq.length() - keylen)
                    return hitCount>0;
                pos++;
                key = Kmer::seq2uint64(seq, pos, keylen-1, valid);
                while(valid == false) {
                    pos++;
                    key = Kmer::seq2uint64(seq, pos, keylen-1, valid);
                    // reach the tail
                    if(pos >= seq.length() - keylen)
                        return hitCount>0;
                }
                continue;
        }
        key = (key << blankBits) >> blankBits;
        bool hit = mKmer->add(key);
        if(hit)
            hitCount++;
    }

    return hitCount>0;
}