#ifndef ACCUMULATE_DENOMINATE_H
#define ACCUMULATE_DENOMINATE_H

#include "parsers.h"

//Keep track on number of ancestrally A,C,G,T bases
//that could be called for a mutation in a given sample

struct DenomCounter{
    uint64_t counts[4];
};

struct DenomCriteria{
    uint32_t DP_min;
    uint32_t DP_max;
    uint32_t mutant_strand_min;
    uint32_t ancestral_in_mutant_max;
    uint32_t mutant_in_ancestral_max;
    double MQ_AD_max;
    double Insert_AD_max;
    double Fisher_std_min;
    double Fisher_map_min;
    
};


typedef std::vector<DenomCounter> DenomCounterVector;

class DenomVisitor : public ReadDataVisitor{
public:
    DenomVisitor(   const RefVector& bam_references, 
                    const SamHeader& header,
                    LocalBamToolsUtils::Fasta& idx_ref,
                    SampleMap& samples, 
                    BamAlignment& ali, 
                    int qual_cut,
                    int mapping_cut,
                    DenomCounterVector &denoms,
                    ModelParams& params,
                    DenomCriteria criteria,
                    double prob_cut);

    virtual ~DenomVisitor() { }

    void Visit(const LocalBamToolsUtils::PileupPosition &pileupData);

    bool MeetsCriteria(SiteStatsSummary stats);

private:
    RefVector m_bam_ref;
    SamHeader m_header;
    BamAlignment& m_ali;
    DenomCounterVector& m_denoms;
    double m_prob_cut;
    const ModelParams &m_params;
    int nsamp;

    SequencingFactory sf;
    MutationMatrix m_mut_paths;
    MutationMatrix m_nomut_paths;
    SampleNames descendant_names;
    DenomCriteria m_criteria;
};
#endif // ACCUMULATE_DENOMINATE_H
