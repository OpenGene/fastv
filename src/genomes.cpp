#include "genomes.h"
#include "util.h"
#include "kmer.h"
#include "editdistance.h"

Genomes::Genomes(string faFile, Options* opt)
{
    mFastaReader = new FastaReader(faFile);
    mOptions = opt;
    mFastaReader->readAll();
    init();
    mMissedCount = 0;
    mHitCount = 0;
}

Genomes::~Genomes()
{
    if(mFastaReader) {
        delete mFastaReader;
        mFastaReader = NULL;
    }

    //cerr << "mMissedCount: " << mMissedCount << endl;
    //cerr << "mHitCount: " << mHitCount << endl;
}

void Genomes::init() {
    initLowComplexityKeys();
    map<string, string> genomes = mFastaReader->contigs();
    map<string, string>::iterator iter;
    mGenomeNum = 0;
    for(iter = genomes.begin(); iter != genomes.end() ; iter++) {
        if(iter->second.size() >= 0xFFFFFF) {
            cerr << "fastv only supports genome size up to 16M, skip " << iter->first << " (" << iter->second.size() << " bp)" << endl;
            break;
        }
        mNames.push_back(iter->first);
        mTotalEditDistance.push_back(0);
        mReads.push_back(0);
        mBases.push_back(0);
        mSequences.push_back(iter->second);
        mCoverage.push_back(vector<uint16>(iter->second.length(), 0));
        mEditDistance.push_back(vector<float>(iter->second.length(), 0));
        mGenomeNum++;
        if(mGenomeNum >= 255) {
            cerr << "fastv only supports up to 255 genomes, other genomes will be skipped." << endl;
            break;
        }
    }

    buildKmerTable();
}

void Genomes::initLowComplexityKeys() {
    int keylen = mOptions->kmerKeyLen;
    const char bases[4] = {'A', 'T', 'C', 'G'};

    // we consider a key with only two positions of different base as low complexity kmer
    for(int i=0; i<4; i++) {
        for(int j=0; j<4; j++) {
            for(int k=0; k<4; k++) {
                char origin = bases[i];
                char diff1 = bases[j];
                char diff2 = bases[k];
                for(int p=0; p<keylen; p++) {
                    for(int q=0; q<keylen; q++) {
                        string seq(keylen, origin);
                        seq[p] = diff1;
                        seq[q] = diff2;
                        bool valid;
                        uint64 key = Kmer::seq2uint64(seq, 0, keylen, valid);
                        mLowComplexityKeys.insert(key);
                    }
                }
            }
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
    // dont add low complexity keys
    if(mLowComplexityKeys.find(key) != mLowComplexityKeys.end())
        return;

    uint32 data = packIdPos(id, pos);
    if(mKmerTable.count(key) == 0)
        mKmerTable[key] = list<uint32>();
    mKmerTable[key].push_back(data);
}

bool Genomes::hasKey(uint64 key) {
    bool hit = mKmerTable.find(key) != mKmerTable.end();
    if(hit)
        mHitCount++;
    else
        mMissedCount++;

    return hit;
}

bool Genomes::align(string& seq) {
    vector<MapResult> results(mGenomeNum);

    int keylen = mOptions->kmerKeyLen;
    int blankBits = 64 - 2*keylen;

    int totalMapped = 0;

    bool valid = true;

    uint32 start = 0;
    uint64 key = Kmer::seq2uint64(seq, start, keylen-1, valid);
    while(valid == false) {
        start++;
        key = Kmer::seq2uint64(seq, start, keylen-1, valid);
        // reach the tail
        if(start >= seq.length() - keylen)
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
                    continue;
                pos++;
                key = Kmer::seq2uint64(seq, pos, keylen-1, valid);
                while(valid == false) {
                    pos++;
                    key = Kmer::seq2uint64(seq, pos, keylen-1, valid);
                    // reach the tail
                    if(pos >= seq.length() - keylen)
                        continue;
                }
                continue;
        }
        key = (key << blankBits) >> blankBits;

        // if the first 10 kmers dont match, then sample it by 10 for speed consideration
        if(pos>10 && pos % 10 !=0)
            continue;

        if(hasKey(key)) {
            list<uint32>& gpList = mKmerTable[key];
            list<uint32>::iterator gpIter;
            for(gpIter = gpList.begin(); gpIter != gpList.end(); gpIter++) {
                // unit32 = 8 bits genome id + 24 bits positions
                uint32 gp = *gpIter;
                uint32 genomeID = 0;
                uint32 genomePos = 0;
                unpackIdPos(gp, genomeID,  genomePos);
                if(results[genomeID].mapped == false) {
                    MapResult r = mapToGenome(seq, pos, mSequences[genomeID], genomePos);

                    if(r.mapped || r.ed < results[genomeID].ed)
                        results[genomeID] = r;
                    if(r.mapped)
                        totalMapped++;
                }
            }

            // no need to process again
            if(totalMapped >= mGenomeNum/2) {
                //cerr << "break at pos: " << pos << endl;
                //break;
            }
        }

    }

    bool mapped = false;
    uint32 minED=0x3FFFFF;
    for(int i=0; i<mGenomeNum; i++) {
        mapped |= results[i].mapped;

        if(results[i].mapped) {
            cover(i, results[i].start, results[i].len);
            mTotalEditDistance[i] += results[i].ed;

            if(minED > results[i].ed)
                minED = results[i].ed;
        }
    }

    return mapped;
}

MapResult Genomes::mapToGenome(string& seq, uint32 seqPos, string& genome, uint32 genomePos) {
    MapResult ret;

    if(genomePos < seqPos)
        return ret;

    uint32 gp = genomePos - seqPos;

    if(genome.length() - genomePos < seq.length())
        return ret;

    uint32 hd = hamming_distance(seq.c_str(), seq.length(), genome.c_str() + gp, seq.length());

    uint32 ed = 0;

    // using hamming distance to accelerate computing edit distance
    if(hd<=2)
        ed = hd;
    else
        ed = edit_distance(seq.c_str(), seq.length(), genome.c_str() + gp, seq.length());

    ret.ed = ed;
    ret.start = gp;
    ret.len = seq.length();
    ret.mapped = ed < 10 && ed < seq.length()/4; // TODO: export to options

    return ret;
}

void Genomes::report() {
    cerr << endl << "Coverage of genomes:" << endl;
    for(int i=0; i<mGenomeNum; i++) {
        cerr << mReads[i] << " reads/" << mBases[i] << " bases/" << mTotalEditDistance[i] << " mismatches: " << mNames[i] << endl;
    }
}

void Genomes::cover(int id, uint32 pos, uint32 len) {
    if(id >= mCoverage.size()) {
        error_exit("WRONG id");
    }

    for(uint32 p=pos; p<pos+len && p<mCoverage[id].size(); p++) {
        mCoverage[id][p]++;
    }

    mReads[id]++;
    mBases[id] += len;
}

uint32 Genomes::packIdPos(uint32 id, uint32 position) {
    uint32 data = 0;
    data |= (id << 24);
    data |= position;
    return data;
}

void Genomes::unpackIdPos(uint32 data, uint32& id, uint32& pos) {
    id = data >> 24;
    pos = data & 0xFFFFFF;
}