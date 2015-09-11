//
// Created by steven on 9/9/15.
//

#include "variant_visitor.h"

VariantVisitor::VariantVisitor(const RefVector &bam_references, LocalBamToolsUtils::Fasta &idx_ref,
                               ostream *out_stream, SampleMap &samples, const ModelParams &p,
                               BamAlignment &ali,
                               int qual_cut, int mapping_cut, double prob_cut) :
        ReadDataVisitor(bam_references, idx_ref, samples, p,
                        ali, qual_cut, mapping_cut),
        m_ostream(out_stream), m_prob_cut(prob_cut) {

    MutationMatrix mt = MutationAccumulation(m_params, true);
    m_mut_paths = MutationAccumulation(m_params, false);
    m_nomut_paths = m_mut_paths - mt;

    qual_cut_char = (char) (qual_cut + 33);
    rg_tag.push_back(ZERO_CHAR);
    rg_tag += "RGZ";
}


void VariantVisitor::Visit(const LocalBamToolsUtils::PileupPosition& pileupData) {
//  if (GatherReadData(pileupData) ){
    if (GatherReadDataV2(pileupData) ){

        double prob = TetMAProbability(m_params, sf, site_data, m_mut_paths, m_nomut_paths);

//        *m_ostream << m_bam_references[pileupData.RefId].RefName << '\t'
//            << pileupData.Position << '\t'
//            << current_base << '\t'
//            << prob << '\t' << endl;

        if(prob >= m_prob_cut){
//            double prob_one = TetMAProbOneMutation(m_params, site_data, m_mut_paths, m_nomut_paths);
            double prob_one = TetMAProbOneMutation(m_params, sf, site_data, m_mut_paths, m_nomut_paths);
            *m_ostream << m_bam_references[pileupData.RefId].RefName << '\t'
                << pileupData.Position << '\t'
                << current_base << '\t'
                << prob << '\t'
                << prob_one << '\t'
                << endl;
        }

    }
}
