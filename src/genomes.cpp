#include "genomes.h"
#include "util.h"
#include "kmer.h"
#include "editdistance.h"
#include <sstream>
#include <memory.h>

// we use 512M memory
const int BLOOM_FILTER_LENGTH = (1<<29);

Genomes::Genomes(string faFile, Options* opt)
{
    mFastaReader = new FastaReader(faFile);
    mOptions = opt;
    mBloomFilterArray = NULL;
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
    if(mBloomFilterArray) {
        delete mBloomFilterArray;
        mBloomFilterArray = NULL;
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
        if(mGenomeNum >= 255) {
            cerr << "fastv only supports up to 255 genomes, other genomes will be skipped." << endl;
            break;
        }
        if(iter->second.size() >= 0xFFFFFF) {
            cerr << "fastv only supports genome size up to 16M, skip " << iter->first << " (" << iter->second.size() << " bp)" << endl;
            continue;
        }
        mNames.push_back(iter->first);
        mTotalEditDistance.push_back(0);
        mReads.push_back(0);
        mBases.push_back(0);
        mSequences.push_back(iter->second);

        mGenomeNum++;
    }

    if(mOptions->statsBinSize == 0) {
        initBinSize();
    }

    for(iter = genomes.begin(); iter != genomes.end() ; iter++) {
        int binNum = (iter->second.length() + 1)/mOptions->statsBinSize;
        mCoverage.push_back(vector<float>(binNum, 0));
        mEditDistance.push_back(vector<float>(binNum, 0));
    }

    buildKmerTable();
    initBloomFilter();
}

void Genomes::initBinSize() {
    int maxSize = 0;
    for(int i=0; i<mGenomeNum; i++) {
        if(mSequences[i].length() > maxSize)
            maxSize = mSequences[i].length();
    }

    int binSize = maxSize / 1600;

    if(binSize < 1)
        mOptions->statsBinSize = 1;
    else if(binSize < 10)
        mOptions->statsBinSize = binSize;
    else if(binSize < 100)
        mOptions->statsBinSize = (binSize/10) * 10;
    else if(binSize < 1000)
        mOptions->statsBinSize = (binSize/100) * 100;
    else if(binSize < 10000)
        mOptions->statsBinSize = (binSize/1000) * 1000;
    else if(binSize < 100000)
        mOptions->statsBinSize = (binSize/10000) * 10000;
    else
        mOptions->statsBinSize = 100000;
}

void Genomes::initBloomFilter() {
    mBloomFilterArray = new char[BLOOM_FILTER_LENGTH];
    memset(mBloomFilterArray, 0, BLOOM_FILTER_LENGTH * sizeof(char));

    //update bloom filter array
    const uint64 bloomFilterFactors[3] = {1713137323, 371371377, 7341234131};

    unordered_map<uint64, list<uint32>>::iterator iter;
    for(iter = mKmerTable.begin(); iter != mKmerTable.end(); iter++) {
        uint64 key = iter->first;
        for(int b=0; b<3; b++) {
            mBloomFilterArray[(bloomFilterFactors[b] * key) & (BLOOM_FILTER_LENGTH-1) ] = 1;
        }
    }
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
                    bool outterBreak = false;
                    key = Kmer::seq2uint64(seq, pos, keylen-1, valid);
                    while(valid == false) {
                        pos++;
                        key = Kmer::seq2uint64(seq, pos, keylen-1, valid);
                        // reach the tail
                        if(pos >= seq.length() - keylen - polyATailLen) {
                            outterBreak = true;
                            break;
                        }
                    }
                    if(outterBreak)
                        break;

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
    // check bloom filter
    const uint64 bloomFilterFactors[3] = {1713137323, 371371377, 7341234131};
    for(int b=0; b<3; b++) {
        if(mBloomFilterArray[(bloomFilterFactors[b] * key) & (BLOOM_FILTER_LENGTH-1)] == 0 )
            return false;
    }

    bool hit = mKmerTable.find(key) != mKmerTable.end();
    if(hit)
        mHitCount++;
    else
        mMissedCount++;

    return hit;
}

bool Genomes::align(string& seq) {
    vector<vector<MapResult>> results(mGenomeNum);

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
                bool outterBreak = false;
                while(valid == false) {
                    pos++;
                    key = Kmer::seq2uint64(seq, pos, keylen-1, valid);
                    // reach the tail
                    if(pos >= seq.length() - keylen){
                        outterBreak = true;
                        break;
                    }
                }
                if(outterBreak)
                    break;

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
                if(results[genomeID].size() == 0) {
                    MapResult r = mapToGenome(seq, pos, mSequences[genomeID], genomePos);

                    if(r.mapped) {
                        totalMapped++;
                        results[genomeID].push_back(r);
                        while(true) {
                            list<uint32>::iterator gpIterNext = gpIter;
                            gpIterNext++;
                            if(gpIter == gpList.end())
                                break;
                            uint32 gpNext = *gpIterNext;
                            uint32 genomeIDNext = 0;
                            uint32 genomePosNext = 0;
                            unpackIdPos(gpNext, genomeIDNext,  genomePosNext);

                            if(genomeIDNext != genomeID) 
                                break;

                            MapResult rNext = mapToGenome(seq, pos, mSequences[genomeID], genomePosNext);
                            if(rNext.mapped) {
                                results[genomeID].push_back(rNext);
                            }
                            gpIter = gpIterNext;
                        }
                    }
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
    for(int i=0; i<mGenomeNum; i++) {
        if(results[i].size() > 0) {
            mapped  = true;
            float frac = 1.0f/results[i].size();

            uint32 minED=0x3FFFFF;
            for(int p=0; p<results[i].size(); p++) {
                cover(i, results[i][p].start, results[i][p].len, results[i][p].ed, frac);
                if(minED > results[i][p].ed)
                    minED = results[i][p].ed;
            }

            mReads[i]++;
            mBases[i] += results[i][0].len;
            mTotalEditDistance[i] += minED;
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
    ret.mapped = ed <= mOptions->edThreshold && ed < seq.length()/4; // TODO: export to options

    return ret;
}

void Genomes::report() {
    cerr << endl << "Coverage of genomes:" << endl;
    for(int i=0; i<mGenomeNum; i++) {
        cerr << mReads[i] << " reads/" << mBases[i] << " bases/" << mTotalEditDistance[i] << " mismatches: " << mNames[i] << endl;
        for(int j=0; j<mCoverage[i].size(); j++) {
            cerr << mCoverage[i][j]/mOptions->statsBinSize << " ";
        }
        cerr << endl;
        for(int j=0; j<mEditDistance[i].size(); j++) {
            cerr << mEditDistance[i][j]/mOptions->statsBinSize << " ";
        }
        cerr << endl;
    }
}

double Genomes::getCoverageRate(int id) {
    if(mCoverage[id].size() == 0)
        return 0.0;

    int covered = 0;

    for(int b=0; b<mCoverage[id].size(); b++) {
        if(mCoverage[id][b] >= mOptions->depthThreshold)
            covered++;
    }

    return (double)covered/(double)mCoverage[id].size();
}

void Genomes::reportJSON(ofstream& ofs) {
    ofs << "\t" << "\"genome_mapping_result\": {" << endl;
    ofs << "\t\t" << "\"genome_number\": " << mGenomeNum << "," << endl;
    ofs << "\t\t" << "\"bin_size\": " << mOptions->statsBinSize << "," << endl;
    ofs << "\t\t" << "\"genome_coverage\": [";
    for(int i=0; i<mGenomeNum; i++) {
        string name = mNames[i];
        long reads = mReads[i];
        long bases = mBases[i];
        long totalED = mTotalEditDistance[i];
        double coverageRate = getCoverageRate(i);

        if(i != 0) 
            ofs << ", " << endl;
        ofs << "\t\t\t{" << endl;
        ofs << "\t\t\t\t\"name\":\"" <<  name << "\"," << endl;
        ofs << "\t\t\t\t\"size\":" <<  mSequences[i].length() << "," << endl;
        ofs << "\t\t\t\t\"reads\":" <<  reads << "," << endl;
        ofs << "\t\t\t\t\"bases\":" <<  bases << "," << endl;
        ofs << "\t\t\t\t\"coverage_rate\":" <<  coverageRate << "," << endl;
        if(bases == 0)
            ofs << "\t\t\t\t\"avg_mismatch_ratio\":" <<  0.0 << "," << endl;
        else
            ofs << "\t\t\t\t\"avg_mismatch_ratio\":" <<  totalED / (double)bases << "," << endl;
        ofs << "\t\t\t\t\"coverage\":[" <<  getCoverageY(i) << "]," << endl;
        ofs << "\t\t\t\t\"mismatch_ratios\":[" <<  getEditDistanceY(i) << "]" << endl;
        ofs << "\t\t\t}";
    }
    ofs << "\t\t]" << endl;
    ofs << "\t}," << endl;
}

void Genomes::reportHtml(ofstream& ofs) {
    ofs << "<div id='genome_coverage' style='display:none;color:white;padding:5px;background-color: rgba(0,0,0,0.6);border:1px dotted #666666;font-size:12px;line-height:15px;'> </div>" << endl;
    ofs << "<script src='http://opengene.org/fastv/coverage.js'></script>" << endl;
    ofs << "<script language='javascript'>" << endl;

    ofs << "var genome_sizes = [";
    for(int i=0; i<mGenomeNum; i++) {
        if(i != 0) 
            ofs << ", ";
        ofs << mSequences[i].length();
    }
    ofs << "];" << endl;

    ofs << "var genome_coverage_data = [" << endl;
    for(int i=0; i<mGenomeNum; i++) {
        string name = mNames[i];
        long reads = mReads[i];
        long bases = mBases[i];
        long totalED = mTotalEditDistance[i];
        double coverageRate = getCoverageRate(i);
        if(i != 0) 
            ofs << ", " << endl;
        ofs << "{" << endl;
        ofs << "\"name\":\"" <<  name << "\"," << endl;
        ofs << "\"reads\":" <<  reads << "," << endl;
        ofs << "\"bases\":" <<  bases << "," << endl;
        ofs << "\"coverage_rate\":" <<  coverageRate * 100 << "," << endl;
        if(bases == 0)
            ofs << "\"avg_mismatch_ratio\":" <<  0.0 << "," << endl;
        else
            ofs << "\"avg_mismatch_ratio\":" <<  totalED / (double)bases << "," << endl;
        ofs << "\"coverage\":[" <<  getCoverageY(i) << "]," << endl;
        ofs << "\"mismatch_ratios\":[" <<  getEditDistanceY(i) << "]" << endl;
        ofs << "}";
    }
    ofs << "];" << endl;

    ofs << "var stats_bin = " << mOptions->statsBinSize << "; " << endl;

    ofs << "genome_coverage_data.sort(sortByCoverageRate);" << endl;

    ofs << "drawCoverages('genome_coverage', genome_coverage_data, genome_sizes, stats_bin);" << endl;

    ofs << "</script>" << endl;

    ofs << "<div id='maptips' style='display:none;color:white;padding:5px;background-color: rgba(0,0,0,0.6);border:1px dotted #666666;font-size:12px;line-height:15px;'> </div>" << endl;
}

string Genomes::getPlotX(int id) {
    stringstream ss;
    for(int x = 0; x < mCoverage[id].size(); x++) {
        if(x > 0) 
            ss << ",";

        ss << x*mOptions->statsBinSize;
    }
    return ss.str();
}

string Genomes::getCoverageY(int id) {
    stringstream ss;
    for(int x = 0; x < mCoverage[id].size(); x++) {
        if(x > 0) 
            ss << ",";

        if(x < mCoverage[id].size() - 1)
            ss  << mCoverage[id][x] / (double)mOptions->statsBinSize ;
        else if(mSequences[id].length() - x * mOptions->statsBinSize == 0)
            ss << "0.0";
        else
            ss  << mCoverage[id][x] / (mSequences[id].length() - x * mOptions->statsBinSize) ;
    }
    return ss.str();
}

string Genomes::getEditDistanceY(int id) {
    stringstream ss;
    for(int x = 0; x < mCoverage[id].size(); x++) {
        if(x > 0) 
            ss << ",";

        if(mCoverage[id][x] > 0)
            ss  << (double)mEditDistance[id][x] / (double)mCoverage[id][x] ;
        else
            ss << "0.0";
    }
    return ss.str();
}

void Genomes::cover(int id, uint32 pos, uint32 len, uint32 ed, float frac) {
    if(id >= mCoverage.size()) {
        error_exit("WRONG id");
    }

    int leftBin = pos / mOptions->statsBinSize;
    int rightBin = (pos+len) / mOptions->statsBinSize;

    if(leftBin == rightBin) {
        if(leftBin < mCoverage[id].size()) {
            mCoverage[id][leftBin] += len * frac;
            mEditDistance[id][leftBin] += ed * frac;
        }
    } else {
        for(int bin = leftBin; bin<rightBin; bin++) {
            int left, right;
            if(bin == leftBin)
                left = pos;
            else
                left = bin * mOptions->statsBinSize;

            if(bin == right)
                right = pos + len;
            else
                right = (bin+1) * mOptions->statsBinSize;

            float proportion = (right - left)/(float)len;

            if(bin < mCoverage[id].size()) {
                mCoverage[id][bin] += (right - left) * frac;
                mEditDistance[id][bin] += ed * proportion * frac;
            }
        }
    }
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