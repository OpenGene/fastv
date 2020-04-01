#ifndef JSON_REPORTER_H
#define JSON_REPORTER_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "options.h"
#include "stats.h"
#include "filterresult.h"
#include <fstream>
#include "virusdetector.h"

using namespace std;

class JsonReporter{
public:
    JsonReporter(Options* opt);
    ~JsonReporter();

    void setDupHist(int* dupHist, double* dupMeanGC, double dupRate);
    void setInsertHist(long* insertHist, int insertSizePeak);
    void report(VirusDetector* vd, FilterResult* result, Stats* preStats1, Stats* postStats1, Stats* preStats2 = NULL, Stats* postStats2 = NULL);

private:
    Options* mOptions;
    int* mDupHist;
    double* mDupMeanGC;
    double mDupRate;
    long* mInsertHist;
    int mInsertSizePeak;
};


#endif