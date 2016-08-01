//
// Created by steven on 9/9/15.
//

#include "variant_visitor.h"

VariantVisitor::VariantVisitor(const RefVector &bam_references, LocalBamToolsUtils::Fasta &idx_ref,
                               ostream *out_stream, SampleMap &samples, SampleNames descendants, 
                               const ModelParams &model_param, BamAlignment &ali,
                               int qual_cut, int mapping_cut, double prob_cut) :
        ReadDataVisitor(idx_ref, samples,
                        qual_cut, mapping_cut),
        m_bam_references(bam_references), m_ostream(out_stream),
        m_params(model_param), m_prob_cut(prob_cut), sf(m_params) {

    MutationMatrix mt = MutationAccumulation(m_params, true);
    m_mut_paths = MutationAccumulation(m_params, false);
    m_nomut_paths = m_mut_paths - mt;
    descendant_names = descendants;

}


void VariantVisitor::Visit(const LocalBamToolsUtils::PileupPosition &pileupData) {
    if (GatherReadData(pileupData)) {

        double prob = TetMAProbability(m_params, sf, site_data, m_mut_paths, m_nomut_paths);

//        *m_ostream << m_bam_references[pileupData.RefId].RefName << '\t'
//            << pileupData.Position << '\t'
//            << current_base << '\t'
//            << prob << '\t' << endl;

        if (prob >= m_prob_cut) {
//            double prob_one = TetMAProbOneMutation(m_params, site_data, m_mut_paths, m_nomut_paths);
            MutationDescription details = DescribeMutant(m_params, sf, site_data, m_mut_paths, m_nomut_paths);
            SiteStatsSummary stats = CalculateStats(pileupData, details.mutant_line, details.mutant_allele_index);

            *m_ostream << m_bam_references[pileupData.RefId].RefName << '\t'
                << pileupData.Position << '\t'
                << pileupData.Position + 1 << '\t'
                << current_base << '\t'
                << descendant_names[ details.mutant_line ] << '\t'
                << details.from_genotype << "->" << details.to_genotype << '\t'
                << prob << '\t'
                << details.line_prob << '\t'
                << details.genotype_prob << '\t'
                << details.lik << '\t'
                << stats.DP << '\t'            
                << stats.NM_F << '\t'
                << stats.NM_R << '\t'
                << stats.N_minor << '\t'
                << stats.NM_WT << '\t'
                << stats.MQ_AD << '\t'
                << stats.Insert_AD << '\t'
                << stats.FisherStrandBias << '\t'
                << stats.FisherPairBias << '\t'
                << endl;
        }

    }
}
