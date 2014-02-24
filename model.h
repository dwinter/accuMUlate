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
};

struct ModelInput{
    string chr;
    uint64_t position;
    uint16_t reference;
    ReadDataVector all_reads;
};

typedef Eigen::Array4d HaploidProbs;
typedef Eigen::Array<double, 16, 1> DiploidProbs;
typedef Eigen::Array<double, 16, 4> MutationMatrix;



double TetMAProbOneMutation(const ModelParams &params, const ModelInput site_data);
double TetMAProbability(const ModelParams &params, const ModelInput site_data);


#endif
