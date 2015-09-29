//
// Created by steven on 8/27/15.
//

#ifndef ACCUMULATE_DATA_STRUCT_H
#define ACCUMULATE_DATA_STRUCT_H

#include <vector>
#include "Eigen/Dense"


union ReadData {
    uint64_t key;
    uint16_t reads[4];
//    ReadData() {  }
//
//    ReadData(uint64_t k) {
//        key = k;
//    }
//
//
//    ReadData(ReadData &&other) : key(other.key) {
////        std::cout << "ReadData move constructor" << std::endl;
//    }
//
//    ReadData &operator=(ReadData &&other) {
////        std::cout << "ReadData move assignment" << std::endl;
//        key = other.key;
//        return *this;
//    }
//
//    ReadData(const ReadData &other) = default;
//
//    ReadData &operator=(const ReadData &other) = default;
};

typedef std::vector<ReadData> ReadDataVector;

struct ModelParams{
    double theta;               //
    std::vector<double> nuc_freq;    //ACGT
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

typedef std::array<double, 10> DiploidProbsIndex10;
typedef std::array<double, 10> Array10D;
typedef std::array<double, 4> Array4D;


#endif //ACCUMULATE_DATA_STRUCT_H
