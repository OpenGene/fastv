#ifndef VIRUSDETECTOR_H
#define VIRUSDETECTOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "options.h"
#include "read.h"
#include "kmer.h"
#include "genomes.h"

using namespace std;

class VirusDetector{
public:
    VirusDetector(Options* opt);
    ~VirusDetector();
    bool detect(Read* r);
    void report();

private:
    Options* mOptions;
    Genomes* mGenomes;
    Kmer* mKmer;
    uint64 mHits;
};


#endif