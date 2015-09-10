#ifndef model_H
#define model_H

#include <src/mutations/sequencing_factory.h>
#include "Eigen/Dense"
#include "data_struct.h"

using namespace std;

TransitionMatrix F81(const ModelParams &params);

MutationMatrix MutationAccumulation(const ModelParams &params, bool and_mut);

/*Version 1*/
DiploidProbs DiploidPopulation(const ModelParams &params, int ref_allele);

DiploidProbs DiploidSequencing(const ModelParams &params, int ref_allele, ReadData data);

HaploidProbs HaploidSequencing(const ModelParams &params, int ref_allele, ReadData data);

GenotypeProbs PopulationProbs(const ModelParams &params, int ref_allele);

GenotypeProbs Sequencing(const ModelParams &params, int ref_allele, ReadData data, int ploidy);

double TetMAProbability(const ModelParams &params, const ModelInput site_data, const MutationMatrix m,
                        const MutationMatrix mn);

double TetMAProbOneMutation(const ModelParams &params, const ModelInput site_data, const MutationMatrix m,
                            const MutationMatrix mn);

/*Version 2*/
GenotypeProbs PopulationProbs(SequencingFactory &sf, int ref_allele, int ploidy_ancestor);

GenotypeProbs Sequencing(SequencingFactory &sf, ReadData data, int ploidy);

double TetMAProbability(const ModelParams &params, SequencingFactory &sf, const ModelInput &site_data,
                        const MutationMatrix &m, const MutationMatrix &mn);

double TetMAProbOneMutation(const ModelParams &params, SequencingFactory &sf, const ModelInput &site_data,
                            const MutationMatrix &m, const MutationMatrix &mn);

#endif
