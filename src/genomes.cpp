#include "genomes.h"
#include "util.h"
#include "kmer.h"

Genomes::Genomes(string faFile, Options* opt)
{
    mFastaReader = new FastaReader(faFile);
    mOptions = opt;
    mFastaReader->readAll();
}

Genomes::~Genomes()
{
    if(mFastaReader) {
        delete mFastaReader;
        mFastaReader = NULL;
    }
}

void Genomes::init() {
    map<string, string> genomes = mFastaReader->contigs();
    map<string, string>::iterator iter;
    int total = 0;
    for(iter = genomes.begin(); iter != genomes.end() ; iter++) {
        mNames.push_back(iter->first);
        mSequences.push_back(iter->second);
        mCoverage.push_back(vector<uint16>(iter->second.length(), 0));
        total++;
        if(total >= 1000) {
            cerr << "fastv only supports up to 1000 genomes, other genomes will be skipped." << endl;
            break;
        }
    }
}

void Genomes::buildKmerTable() {
    int keylen = mOptions->kmerKeyLen;
    int blankBits = 64 - 2*keylen;
    const int polyATailLen = 28;
    bool valid = true;
    for(uint32 i=0; i<mNames.size(); i++) {
        string& seq = mSequences[i];
        if(seq.length() < keylen)
            continue;
        // first calculate the first keylen-1 kmer
        // skip the polyA tail
        uint32 start = 0;
        uint64 key = Kmer::seq2uint64(seq, start, keylen-1, valid);
        while(valid == false) {
            start++;
            key = Kmer::seq2uint64(seq, start, keylen-1, valid);
            // reach the tail
            if(start >= seq.length() - keylen - polyATailLen)
                return;
        }
        for(uint32 pos = start; pos < seq.length() - keylen - polyATailLen; pos++) {
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
                    pos++;
                    key = Kmer::seq2uint64(seq, pos, keylen-1, valid);
                    while(valid == false) {
                        pos++;
                        key = Kmer::seq2uint64(seq, pos, keylen-1, valid);
                        // reach the tail
                        if(pos >= seq.length() - keylen - polyATailLen)
                            break;
                    }
                    continue;
            }
            key = (key << blankBits) >> blankBits;
            addKmer(key, i, pos);
        }
    }
}

void Genomes::addKmer(uint64 key, uint32 id, uint32 pos) {
    uint32 data = packIdPos(id, pos);
    if(mKmerTable.count(key) == 0)
        mKmerTable[key] = list<uint32>();
    mKmerTable[key].push_back(data);
}

void Genomes::cover(int id, uint32 pos, uint32 len) {
    if(id >= mCoverage.size()) {
        error_exit("WRONG id");
    }

    for(uint32 p=pos; p<pos+len && p<mCoverage[id].size(); p++) {
        mCoverage[id][p]++;
    }
}

uint32 Genomes::packIdPos(uint32 id, uint32 position) {
    uint32 data = 0;
    data |= (id << 22);
    data |= position;
    return data;
}

void Genomes::unpackIdPos(uint32 data, uint32& id, uint32& pos) {
    id = data >> 22;
    pos = data & 0x3FFFFF;
}