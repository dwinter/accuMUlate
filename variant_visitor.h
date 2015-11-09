//
// Created by steven on 9/9/15.
//

#ifndef ACCUMULATE_VARIANT_VISITOR_H
#define ACCUMULATE_VARIANT_VISITOR_H

#include "parsers.h"

class VariantVisitor : public ReadDataVisitor {
public:

    VariantVisitor(const RefVector &bam_references,
                   LocalBamToolsUtils::Fasta &idx_ref,
                   ostream *out_stream,
                   SampleMap &samples,
                   const ModelParams &model_param,
                   BamAlignment &ali,
                   int qual_cut, int mapping_cut, double prob_cut);

    virtual ~VariantVisitor() { }

    void Visit(const LocalBamToolsUtils::PileupPosition &pileupData);

private:
    const RefVector &m_bam_references;
    ostream *m_ostream;
    const ModelParams &m_params;
    double m_prob_cut;

    SequencingFactory sf;
    MutationMatrix m_mut_paths;
    MutationMatrix m_nomut_paths;

};


#endif //ACCUMULATE_VARIANT_VISITOR_H
