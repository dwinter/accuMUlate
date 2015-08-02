#ifndef parsers_H
#define parsers_H

#include "utils/bamtools_pileup_engine.h"
#include "model.h"
#include <unordered_map>

using namespace std;
using namespace BamTools;

//typedef vector< string > SampleNames;
typedef unordered_map<string, uint32_t> SampleMap;

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
        ReadDataVisitor(const RefVector& bam_references, 
                        BamTools::Fasta& idx_ref,
                        SampleMap& samples, 
                        const ModelParams& p,  
                        BamAlignment& ali, 
                        int qual_cut,
                        int mapping_cut);
    public: 
        bool GatherReadData(const PileupPosition& pileupData) ;
    public:
        char current_base;
        string tag_id;
        uint64_t chr_index;
        ModelInput site_data;
        const ModelParams& m_params;
        const RefVector& m_bam_references;
        MutationMatrix m_mutation_paths;
        MutationMatrix m_non_mutation_paths;

    private:
        //set by arguments
        BamTools::Fasta& m_idx_ref;
        SampleMap& m_samples;
        BamAlignment& m_ali;
        int m_qual_cut;
        int m_mapping_cut;
        //refered to by fnxs
};
        
 
//ModelInput CollectReadData(BamTools::PileupPosition& pileupData);



//Helper functions

bool include_site(BamTools::PileupAlignment pileup, uint16_t map_cut, uint16_t qual_cut);
uint16_t base_index(char b);
string get_sample(string& tag);
//uint32_t find_sample_index(string, SampleNames);

#endif



