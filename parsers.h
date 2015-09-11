#ifndef parsers_H
#define parsers_H

#include "src/io_data/local_bamtools/bamtools_pileup_engine.h"
#include "src/io_data/local_bamtools/bamtools_fasta.h"
//#include "utils/bamtools_pileup_engine.h"
//#include "utils/bamtools_fasta.h"

#include "model.h"
#include <unordered_map>
#include <src/mutations/sequencing_factory.h>

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

class ReadDataVisitor : public LocalBamToolsUtils::PileupVisitor{

    public:
        static const char ZERO_CHAR = ((char) 0);
        ReadDataVisitor(const RefVector& bam_references,
                        LocalBamToolsUtils::Fasta& idx_ref,
                        SampleMap& samples,
                        const ModelParams& p,
                        BamAlignment& ali,
                        int qual_cut,
                        int mapping_cut);

    virtual ~ReadDataVisitor() { }

    public:

        [[deprecated]]
        bool GatherReadData(const LocalBamToolsUtils::PileupPosition& pileupData) ;

        bool GatherReadDataV2(const LocalBamToolsUtils::PileupPosition &pileupData) ;
        uint32_t GetSampleIndex(const string &tag_data);

    protected:
        const ModelParams& m_params;
        const RefVector& m_bam_references;

        std::string rg_tag;
        ModelInput site_data;
        SequencingFactory sf;
        char qual_cut_char;
        char current_base;

        uint32_t total_sample_count;

    //        uint64_t chr_index;
        [[deprecated]]
        string tag_id;
        [[deprecated]]
        MutationMatrix m_mutation_paths;
        [[deprecated]]
        MutationMatrix m_non_mutation_paths;

    private:
        //set by arguments
        LocalBamToolsUtils::Fasta& m_idx_ref;
        SampleMap& m_samples;
        BamAlignment& m_ali;
        int m_qual_cut;
        int m_mapping_cut;
        //refered to by fnxs




};
        
 
//ModelInput CollectReadData(BamTools::PileupPosition& pileupData);

//Helper functions
//string get_sample(string& tag);
//uint32_t find_sample_index(string, SampleNames);
[[deprecated]]
uint16_t base_index(char b);
[[deprecated]]
bool include_site(LocalBamToolsUtils::PileupAlignment pileup, uint16_t map_cut, uint16_t qual_cut);

//Update version
bool include_site_v2(const BamTools::BamAlignment &alignment, const int &pos, const uint16_t &map_cut,
                     const char &qual_cut);
extern const int base_index_lookup[128];



#endif



