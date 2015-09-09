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


    DiploidProbs GetDiploidSequencing(int ref_allele, ReadData &data);

    HaploidProbs GetHaploidSequencing(int ref_allele, ReadData &data);

private:

    ModelParams model_params;
    double phi_haploid;
    double phi_diploid;
    double error_prob;
    double theta;
    std::vector<double> frequency_prior;
    Array10D ancestor_prior;

    double haploid_alphas[4][4];

    std::vector<SiteGenotypesIndex> sgi;

private:
    std::array<DiploidProbs, 4> ref_diploid_probs;
    HaploidProbs ref_halpoid_probs;


    std::vector<HaploidProbs> convert_index_key_to_haploid;
    std::vector<DiploidProbsIndex10> convert_index_key_to_diploid_10;

    std::vector<double> convert_index_key_to_haploid_scaler;
    std::vector<double> convert_index_key_to_diploid_10_scaler;

//    std::vector<HaploidProbs> convert_index_key_to_haploid_unnormalised;
//    std::vector<DiploidProbsIndex10> convert_index_key_to_diploid_10_unnormalised;


    std::unordered_map<uint64_t, uint32_t> map_rd_to_index;
    std::array<std::unordered_map<uint64_t, uint32_t>, 4> map_ancestor_to_index;

    std::unordered_map<int, int> map_des_count;

    //
//    void CalculateDescendantGenotypes(SiteGenotypes &seq_prob);
//    void CalculateAncestorGenotype(SiteGenotypes &seq_prob);
//
//    void CalculateDescendantGenotypesIndex(SiteGenotypesIndex &seq_prob);
//    void CalculateAncestorGenotypeIndex(SiteGenotypesIndex &seq_prob);
    void CalculateAncestorPrior();

    DiploidProbs CreateRefDiploidProbs(int ref_allele);

    DiploidProbs DiploidSequencing(ReadData const &data);

    HaploidProbs HaploidSequencing(ReadData const &data);

    DiploidProbsIndex10 ConvertDiploid16ToDiploid10(DiploidProbs probs, int reference);


    DiploidProbs CalculateGenotypeProbCache(int ref_allele, ReadData data, int ploidy);


    int index_descendant;
};
#endif //_ACCUMULATE_SEQUENCING_FACTORY_H_
