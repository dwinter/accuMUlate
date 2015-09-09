//
// Created by steven on 3/31/15.
//
#pragma once
#ifndef _SEQUENCING_FACTORY_H_
#define _SEQUENCING_FACTORY_H_

//#include "model.h"
#include "../lookup.h" //HACK: FIXME: do it properly later
#include "site_genotypes.h"
#include "site_genotype_index.h"
#include "data_struct.h"
#include <unordered_map>
#include "src/distributions/DirichletMultinomialDistribution.h"

class SequencingFactory {


public:

    [[deprecated]]
    typedef std::vector<ModelInput> GenomeData;

    SequencingFactory() = default;

    SequencingFactory(ModelParams const &model_params);

    void CreateSequenceProbsVector(std::vector<SiteGenotypesIndex> &sgi, GenomeData &data);

    void CreateSequenceProbsVector(std::vector<SiteGenotypesIndex> &sgi, ModelInput &data);

    void CreateSequenceProbsVector(GenomeData &data);

    const std::vector<HaploidProbs> RemoveConvertIndexKeyToHaploid();

    const std::vector<DiploidProbsIndex10> RemoveConvertIndexKeyToDiploidIndex10();

    std::vector<double> &&RemoveConvertIndexKeyToDiploidIndex10Scaler();

    std::vector<double> &&RemoveConvertIndexKeyToHaploidScaler();

    std::vector<SiteGenotypesIndex> &&RemoveSiteGenotypeIndexVector();

//    const std::vector<HaploidProbs> RemoveConvertIndexKeyToHaploidUnnormalised();
//    const std::vector<DiploidProbsIndex10> RemoveConvertIndexKeyToDiploidIndex10Unnormalised();

//    void CreateSequenceProbsVector(std::vector<SiteGenotypes> &sp, GenomeData &data);

//    void CreateSequenceProbV1(std::vector<SequenceProb> &sp, GenomeData &data);


//    std::vector<DiploidProbs> &GetConvertIndexKeyToDiploid();
//    std::array<DiploidProbs, 4> &GetRefDiploidProbs();

public:
    const DiploidProbs &getRefDiploidProbs(int ref) const {
        return ref_diploid_probs[ref];
    }

    const HaploidProbs &getRefHaploidProbs() const {
        return ref_halpoid_probs;
    }


    DiploidProbs GetDiploidSequencing(ReadData &data);

    HaploidProbs GetHaploidSequencing(ReadData &data);

private:

    ModelParams model_params;
    double phi_haploid;
    double phi_diploid;
    double error_prob;
    double theta;
    std::vector<double> frequency_prior;


    double haploid_alphas[4][4];
    double alphas_total_haploid;
    double alphas_total_diploid;

    [[deprecated]] Array10D ancestor_prior;

private:
    std::array<DiploidProbs, 4> ref_diploid_probs;
    HaploidProbs ref_halpoid_probs;

    std::vector<HaploidProbs> convert_index_key_to_haploid;
    std::vector<DiploidProbs> convert_index_key_to_diploid;

    std::unordered_map<uint64_t, uint32_t> map_rd_to_index;
    std::unordered_map<uint64_t, uint32_t> map_ancestor_to_index;

    [[deprecated]] std::vector<DiploidProbsIndex10> convert_index_key_to_diploid_10;
    [[deprecated]] std::vector<double> convert_index_key_to_haploid_scaler;
    [[deprecated]] std::vector<double> convert_index_key_to_diploid_10_scaler;
    [[deprecated]] std::array<std::unordered_map<uint64_t, uint32_t>, 4> map_ancestor_ref_to_index_;
    [[deprecated]] std::unordered_map<int, int> map_des_count;
    [[deprecated]] std::vector<SiteGenotypesIndex> sgi;

    //
//    void CalculateDescendantGenotypes(SiteGenotypes &seq_prob);
//    void CalculateAncestorGenotype(SiteGenotypes &seq_prob);
//
//    void CalculateDescendantGenotypesIndex(SiteGenotypesIndex &seq_prob);
//    void CalculateAncestorGenotypeIndex(SiteGenotypesIndex &seq_prob);

    [[deprecated]]
    void CalculateAncestorPrior();
    [[deprecated]]
    DiploidProbsIndex10 ConvertDiploid16ToDiploid10(DiploidProbs probs, int reference);


    DiploidProbs CreateRefDiploidProbs(int ref_allele);

    DiploidProbs DiploidSequencing(ReadData const &data);

    HaploidProbs HaploidSequencing(ReadData const &data);


    int index_descendant;
    int index_ancestor;

public:
    int count1 = 0;
    int count2 = 0;

};
#endif //_ACCUMULATE_SEQUENCING_FACTORY_H_
