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

//    SequencingFactory() = default;
    SequencingFactory(ModelParams const &model_params);

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


    int index_descendant;
    int index_ancestor;
    double haploid_alphas[4][4];
    double alphas_total_haploid;
    double alphas_total_diploid;

    std::array<DiploidProbs, 4> ref_diploid_probs;
    HaploidProbs ref_halpoid_probs;

    std::vector<HaploidProbs> convert_index_key_to_haploid;
    std::vector<DiploidProbs> convert_index_key_to_diploid;

    std::unordered_map<uint64_t, uint32_t> map_rd_to_index;
    std::unordered_map<uint64_t, uint32_t> map_ancestor_to_index;


    DiploidProbs CreateRefDiploidProbs(int ref_allele);

    DiploidProbs DiploidSequencing(ReadData const &data);

    HaploidProbs HaploidSequencing(ReadData const &data);


};
#endif //_ACCUMULATE_SEQUENCING_FACTORY_H_
