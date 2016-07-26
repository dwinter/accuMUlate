#ifndef parsers_H
#define parsers_H

#include <unordered_map>

#include "src/io_data/local_bamtools/bamtools_pileup_engine.h"
#include "src/io_data/local_bamtools/bamtools_fasta.h"
#include "model.h"


using namespace std;
using namespace BamTools;

struct BedInterval {
    string chr;
    uint64_t start;
    uint64_t end;
};

class BedFile {
//        string bed_file_name;
public:
    BedFile(string bed_file_name);
    int get_interval(BedInterval &current_interval);
    ifstream bed_file;

};


class ReadDataVisitor : public LocalBamToolsUtils::PileupVisitor {

public:
    const char ZERO_CHAR = ((char) 0);
    const std::string RG_TAG{{ZERO_CHAR, 'R', 'G', 'Z'}};
    uint64_t region_start = 0;
    uint64_t region_end = -1;
    ReadDataVisitor(LocalBamToolsUtils::Fasta &idx_ref, SampleMap &samples,
                    int qual_cut, int mapping_cut);

    virtual ~ReadDataVisitor() { }

    bool GatherReadData(const LocalBamToolsUtils::PileupPosition &pileupData);

    void SetRegion( BedInterval target_region );

    uint32_t GetSampleIndex(const string &tag_data);


protected:
    ModelInput site_data;
    char current_base;

private:
    LocalBamToolsUtils::Fasta &m_idx_ref;
    SampleMap &m_samples;
    int m_mapping_cut;
    char qual_cut_char;
    uint32_t total_sample_count;

};


//Helper functions
//Update version
bool include_site(const BamTools::BamAlignment &alignment, const int &pos, const uint16_t &map_cut,
                  const char &qual_cut);

extern const int base_index_lookup[128];

//ModelInput CollectReadData(BamTools::PileupPosition& pileupData);
//string get_sample(string& tag);
//uint32_t find_sample_index(string, SampleNames);



#endif



