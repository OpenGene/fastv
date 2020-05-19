
#include "fastareader.h"
#include "util.h"
#include <sstream>
#include <string.h>

FastaReader::FastaReader(string faFile, bool forceUpperCase)
{
    // Set locale and disable stdio synchronization to improve iostream performance
    // http://www.drdobbs.com/the-standard-librarian-iostreams-and-std/184401305
    // http://stackoverflow.com/questions/5166263/how-to-get-iostream-to-perform-better
    setlocale(LC_ALL,"C");
    ios_base::sync_with_stdio(false);

    mFilename = faFile;
    mForceUpperCase = forceUpperCase;

    if (ends_with(mFilename, ".fasta.gz") || ends_with(mFilename, ".fa.gz") || ends_with(mFilename, ".fna.gz")){
        mZipFile = gzopen(mFilename.c_str(), "r");
        mZipped = true;
    }
    else if(ends_with(mFilename, ".fasta") || ends_with(mFilename, ".fa") || ends_with(mFilename, ".fna")){
        mFile.open(mFilename.c_str(), ifstream::in);
        mZipped = false;
    } else {
        error_exit("FASTA file should have a name (*.fasta, *.fa or *.fna) or (*.fasta.gz, *.fa.gz or *.fna.gz). Not a FASTA file: " + mFilename);
    }

    char c;
    // seek to first contig
    while (getChar(c) && c != '>') {
        if (eof()) {
            break;
        }
    }
}

FastaReader::~FastaReader()
{
    if (mZipped){
        if (mZipFile){
            gzclose(mZipFile);
            mZipFile = NULL;
        }
    }
    else {
        if (mFile.is_open()){
            mFile.close();
        }
    }
}

bool FastaReader::getLine(char* line, int maxLine){
    bool status = true;
    if(mZipped)
        status = gzgets(mZipFile, line, maxLine);
    else {
        mFile.getline(line, maxLine, '\n');
        status = !mFile.fail();
    }

    // trim \n, \r or \r\n in the tail
    int readed = strlen(line);
    if(readed >=2 ){
        if(line[readed-1] == '\n' || line[readed-1] == '\r'){
            line[readed-1] = '\0';
            if(line[readed-2] == '\r')
                line[readed-2] = '\0';
        }
    }

    return status;
}

bool FastaReader::eof() {
    if (mZipped) {
        return gzeof(mZipFile);
    } else {
        return mFile.eof();
    }
}

bool FastaReader::getChar(char& c) {
    bool status = true;
    if (mZipped) {
        c = (char)gzgetc(mZipFile);
        if(c == -1)
            status = false;
    } else {
        mFile.get(c);
        status = !mFile.fail();
    }
    return status;
}

void FastaReader::readNext()
{
    const int maxLine = 1024;
    char linebuf[maxLine];

    mCurrentID = "";
    mCurrentDescription = "";
    mCurrentSequence = "";
    bool foundHeader = false;
    
    char c;
    stringstream ssSeq;
    stringstream ssHeader;
    while(true){
        getChar(c);
        if(c == '>' || eof())
            break;
        else {
            if (foundHeader){
                if(mForceUpperCase && c>='a' && c<='z') {
                    c -= ('a' - 'A');
                }
                ssSeq << c;
            }
            else
                ssHeader << c;
        }
        string line;
        if(mZipped) {
            getLine(linebuf, maxLine);
            line = string(linebuf);
        } else {
            getline(mFile,line,'\n');
        }

        if(foundHeader == false) {
            ssHeader << line;
            foundHeader = true;
        }
        else {
            str_keep_valid_sequence(line, mForceUpperCase);
            ssSeq << line;
        }
    }
    mCurrentSequence = ssSeq.str();
    string header = ssHeader.str();

    mCurrentID = header;
}

bool FastaReader::hasNext() {
    return !eof();
}

void FastaReader::readAll() {
    while(!eof()){
        readNext();
        mAllContigs[mCurrentID] = mCurrentSequence;
    }
}

bool FastaReader::test(){
    FastaReader reader("testdata/tinyref.fa");
    reader.readAll();

    string contig1 = "GATCACAGGTCTATCACCCTATTAATTGGTATTTTCGTCTGGGGGGTGTGGAGCCGGAGCACCCTATGTCGCAGT";
    string contig2 = "GTCTGCACAGCCGCTTTCCACACAGAACCCCCCCCTCCCCCCGCTTCTGGCAAACCCCAAAAACAAAGAACCCTA";

    if(reader.mAllContigs.count("contig1") == 0 || reader.mAllContigs.count("contig2") == 0 )
        return false;

    if(reader.mAllContigs["contig1"] != contig1 || reader.mAllContigs["contig2"] != contig2 )
        return false;

    return true;

}



