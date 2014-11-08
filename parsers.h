#ifndef parsers_H
#define parsers_H

#include "utils/bamtools_pileup_engine.h"
#include <unordered_map>

using namespace std;

//typedef vector< string > SampleNames;
typedef unordered_map<string, uint16_t> SampleMap;

struct FastaReferenceData{
    string name;
    uint32_t length; 
    uint64_t end; //endpoint relative to entire length of whole ref
};

struct BedInterval{
    string chr;
    uint64_t start;
    uint64_t end;
};

typedef vector<FastaReferenceData> FastaReferenceVector;

class FastaReference{
        //string ref_file_name;
    public:
        FastaReference(string ref_file_name);
        FastaReferenceVector chromosomes;
        void get_ref_id(string name, int& chr_id);
    private:
        ifstream ref_file;
};

class BedFile{
//        string bed_file_name;
    public:
        BedFile(string bed_file_name);
        int get_interval(BedInterval& current_interval);
        ifstream bed_file;
        
};
//Helper functions

bool include_site(BamTools::PileupAlignment pileup, uint16_t map_cut, uint16_t qual_cut);
uint16_t base_index(char b);
string get_sample(string& tag);
//uint32_t find_sample_index(string, SampleNames);

#endif



