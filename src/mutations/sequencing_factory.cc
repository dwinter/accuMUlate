//
// Created by steven on 3/31/15.
//


#include "sequencing_factory.h"


SequencingFactory::SequencingFactory(ModelParams const &model_params2) :model_params(model_params2) {

    phi_haploid = model_params.phi_haploid;
    phi_diploid = model_params.phi_diploid;
    error_prob = model_params.error_prob;
    theta = model_params.theta;
	frequency_prior = model_params.nuc_freq;

    index_descendant = 0;
    index_ancestor = 0;

    alphas_total_diploid = (1.0-phi_diploid)/phi_diploid;
    alphas_total_haploid = (1.0-phi_haploid)/phi_haploid;

	for (int i : { 0, 1, 2, 3 }) {
        std::fill(haploid_alphas[i], haploid_alphas[i] + 4, error_prob / 3.0 * alphas_total_haploid);
        for (int k : {0, 1, 2, 3}) {
            if (k == i) {
                haploid_alphas[i][k] = (1.0 - error_prob) * alphas_total_haploid;
                break;
            }
        }
    }

    for( int i :{0,1,2,3}) {
        ref_halpoid_probs[i] = model_params.nuc_freq[i];
    }
    for (int i : { 0, 1, 2, 3 }) {
        ref_diploid_probs[i] = CreateRefDiploidProbs(i);
    }

    CalculateAncestorPrior();
//    cout << ref_diploid_probs << endl;

//    index = 0;
//    index_ancestor = 0;

}


DiploidProbs SequencingFactory::CreateRefDiploidProbs(int ref_allele) {
    ReadData d;
    DiploidProbs result;

    double alphas[4];
    for(int i : {0,1,2,3}){
        alphas[i] = theta * frequency_prior[i];
    }

    for(int i : {0,1,2,3}) {
        for(int j=0;j<i;++j) {
            d.key = 0;
            d.reads[ref_allele] = 1;
            d.reads[i] += 1;
            d.reads[j] += 1;
            result[i+j*4] = DirichletMultinomialLogProbability(alphas, d);
            result[j+i*4] = result[i+j*4];
        }
        d.key = 0;
        d.reads[ref_allele] = 1;
        d.reads[i] += 2;
        result[i+i*4] = DirichletMultinomialLogProbability(alphas, d);
    }
    return result.exp();
}


DiploidProbs SequencingFactory::DiploidSequencing(ReadData const &data) {
    DiploidProbs result;
    //TODO: Refactor these alphas
    for(int i : {0,1,2,3}) {
        for(int j=0;j<i;++j) {
            double alphas[4];
            for(int k : {0,1,2,3}) {
                if(k == i || k == j)
                    alphas[k] = (0.5-error_prob/3.0)*alphas_total_diploid;
                else
                    alphas[k] = (error_prob/3.0)*alphas_total_diploid;
            }
            result[i*4+j] = DirichletMultinomialLogProbability(alphas, data);
            result[j*4+i] = result[i*4+j];
        }
        double alphas[4];
        for(int k : {0,1,2,3}) {
            if(k == i)
                alphas[k] = (1.0-error_prob)*alphas_total_diploid;
            else
                alphas[k] = error_prob/3.0*alphas_total_diploid;
        }
        result[i*4+i] = DirichletMultinomialLogProbability(alphas, data);
    }
    double scale = result.maxCoeff();
    return (result - scale).exp();
}



HaploidProbs SequencingFactory::HaploidSequencing(ReadData const &data) {
    HaploidProbs result;
    for (int i : { 0, 1, 2, 3 }) {
        result[i] = DirichletMultinomialLogProbability(haploid_alphas[i], data);
    }
    double scale = result.maxCoeff();
    return (result - scale).exp();
}


void SequencingFactory::CalculateAncestorPrior() {
    for (int i = 0; i < 4; ++i) {
        for (int j = i; j < 4; ++j) {
            int index10 = LookupTable::index_converter_4_4_to_10[i][j];
            ancestor_prior[index10] = frequency_prior[i] * frequency_prior[j];
            if(i != j){
                ancestor_prior[index10] *= 2; //Count both AC and CA
            }
        }
    }
}

const std::vector<HaploidProbs> SequencingFactory::RemoveConvertIndexKeyToHaploid() {
    return std::move(convert_index_key_to_haploid);
}

const std::vector<DiploidProbsIndex10> SequencingFactory::RemoveConvertIndexKeyToDiploidIndex10() {
    return std::move(convert_index_key_to_diploid_10);
}





DiploidProbs SequencingFactory::GetDiploidSequencing(ReadData &data) {

    auto rd_key = data.key;
    auto find_key = map_ancestor_to_index.find(rd_key);
    if (find_key == map_ancestor_to_index.end()) {

        DiploidProbs diploid_probs = DiploidSequencing(data);
        convert_index_key_to_diploid.push_back(diploid_probs);

        map_ancestor_to_index[rd_key] = index_ancestor;
        index_ancestor++;
    }
    auto prob = convert_index_key_to_diploid[map_ancestor_to_index[rd_key]];
    return prob;
}

HaploidProbs SequencingFactory::GetHaploidSequencing(ReadData &data) {

    auto rd_key = data.key;
    auto find_key = map_rd_to_index.find(rd_key);
    if (find_key == map_rd_to_index.end()) {

        HaploidProbs haploid_prob = HaploidSequencing(data);
        convert_index_key_to_haploid.push_back(haploid_prob);

        map_rd_to_index[rd_key] = index_descendant;
        index_descendant++;
    }
    auto prob = convert_index_key_to_haploid[map_rd_to_index[rd_key]];
    return prob;
}
