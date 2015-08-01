#ifndef model_H
#define model_H

#include "Eigen/Dense"

using namespace std;

union ReadData{
    uint16_t reads[4];
    uint64_t key;
};

typedef vector<ReadData> ReadDataVector;

struct ModelParams{
    double theta;               //
    vector<double> nuc_freq;    //ACGT
    double mutation_rate;       //
    double error_prob;          // Sequencing error-rate 
    double phi_haploid;         // Overdispersion for haploid sequencing
    double phi_diploid;         // Overdispersion for diploid sequencing
    int ploidy_ancestor;         //
    int ploidy_descendant;       //
};

struct ModelInput{// Can probably stand to lose this, started out more complex..
    uint16_t reference;
    ReadDataVector all_reads;
};

typedef Eigen::ArrayXd GenotypeProbs;

typedef Eigen::Array4d HaploidProbs;
typedef Eigen::Array<double, 16, 1> DiploidProbs;

typedef Eigen::Matrix4d TransitionMatrix;

typedef Eigen::ArrayXXd MutationMatrix; //array of dynamic dimensions

//typedef Eigen::Array<double, 4, 4> MutationMatrix_HH;
//typedef Eigen::Array<double, 16, 4> MutationMatrix_DH;
//typedef Eigen::Array<double, 16, 16> MutationMatrix_DD;


DiploidProbs DiploidSequencing(const ModelParams &params, int ref_allele, ReadData data); 
double TetMAProbability(const ModelParams &params, const ModelInput site_data, const MutationMatrix m, const MutationMatrix mn);
double TetMAProbOneMutation(const ModelParams &params, const ModelInput site_data, const MutationMatrix m, const MutationMatrix mn);
MutationMatrix MutationAccumulation(const ModelParams &params, bool and_mut);

#endif
