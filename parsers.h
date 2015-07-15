#ifndef parsers_H
#define parsers_H

#include "utils/bamtools_pileup_engine.h"
#include "model.h"
#include <unordered_map>

using namespace std;

//typedef vector< string > SampleNames;
typedef unordered_map<string, uint16_t> SampleMap;

struct BedInterval{
    string chr;
    uint64_t start;
    uint64_t end;
};

class BedFile{
//        string bed_file_name;
    public:
        BedFile(string bed_file_name);
        int get_interval(BedInterval& current_interval);
        ifstream bed_file;
        
};

class ReadDataVisitor : public BamTools::PileupVisitor{
    public: 
        bool GatherReadData() {};
};
        
 
//ModelInput CollectReadData(BamTools::PileupPosition& pileupData);



//Helper functions

bool include_site(BamTools::PileupAlignment pileup, uint16_t map_cut, uint16_t qual_cut);
uint16_t base_index(char b);
string get_sample(string& tag);
//uint32_t find_sample_index(string, SampleNames);

#endif



