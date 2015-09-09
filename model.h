#ifndef model_H
#define model_H

#include <src/mutations/sequencing_factory.h>
#include "Eigen/Dense"
#include "data_struct.h"

using namespace std;

//typedef Eigen::Array<double, 4, 4> MutationMatrix_HH;
//typedef Eigen::Array<double, 16, 4> MutationMatrix_DH;
//typedef Eigen::Array<double, 16, 16> MutationMatrix_DD;
DiploidProbs DiploidSequencing(const ModelParams &params, int ref_allele, ReadData data);

double TetMAProbability(const ModelParams &params, const ModelInput site_data, const MutationMatrix m,
                        const MutationMatrix mn);

double TetMAProbOneMutation(const ModelParams &params, const ModelInput site_data, const MutationMatrix m,
                            const MutationMatrix mn);

MutationMatrix MutationAccumulation(const ModelParams &params, bool and_mut);

//NEW
double TetMAProbability(const ModelParams &params, SequencingFactory &sf, const ModelInput &site_data,
                        const MutationMatrix &m, const MutationMatrix &mn);

double TetMAProbOneMutation(const ModelParams &params, SequencingFactory &sf, const ModelInput site_data,
                            const MutationMatrix m, const MutationMatrix mn);

GenotypeProbs PopulationProbs(SequencingFactory &sf, int ref_allele, int ploidy_ancestor);

GenotypeProbs Sequencing(SequencingFactory &sf, ReadData data, int ploidy);


#endif
