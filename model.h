#ifndef model_H
#define model_H

#include "src/mutations/sequencing_factory.h"
#include "data_struct.h"

using namespace std;

TransitionMatrix F81(const ModelParams &params);

MutationMatrix MutationAccumulation(const ModelParams &params, bool and_mut);

GenotypeProbs PopulationProbs(SequencingFactory &sf, int ref_allele, int ploidy_ancestor);

GenotypeProbs Sequencing(SequencingFactory &sf, ReadData data, int ploidy);

double TetMAProbability(const ModelParams &params, SequencingFactory &sf, const ModelInput &site_data,
                        const MutationMatrix &m, const MutationMatrix &mn);

double TetMAProbOneMutation(const ModelParams &params, SequencingFactory &sf, const ModelInput &site_data,
                            const MutationMatrix &m, const MutationMatrix &mn);




#endif
