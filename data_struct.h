//
// Created by steven on 8/27/15.
//

#ifndef ACCUMULATE_DATA_STRUCT_H
#define ACCUMULATE_DATA_STRUCT_H

#include <vector>
#include <unordered_map>
#include <iostream>

#include <Eigen/Dense>


const uint16_t MAX_UINT16 = std::numeric_limits<uint16_t>::max();
const uint32_t MAX_UINT32 = std::numeric_limits<uint32_t>::max();

union ReadData {
    uint64_t key;
    uint16_t reads[4];

    ReadData() {
        key=0;
    }

    explicit ReadData(uint64_t k) { //for ReadData (1000)
        key = k;
    }

    ReadData(uint16_t (&r4)[4]) { // for ReadData (array_var)
        std::copy(std::begin(r4), std::end(r4), std::begin(reads));
//        std::memcpy(reads, r4, sizeof(reads));
    }

    ReadData(std::initializer_list<uint64_t > r4){
        if(r4.size() == 1){ //for ReadData {10000000}
            key = *r4.begin();
        }
        else if (r4.size() == 4){ //for ReadData{{1,2,3,4}}
            int i = 0;
            for (auto item : r4) {
                if(item > MAX_UINT16){
                    std::cerr << "ReadData Warning!! " << item << " > " << MAX_UINT16 << std::endl;
                }

                reads[i++] = item;
            }
        }
    }

    ReadData(ReadData &&other) : key(other.key) {
//        std::cout << "ReadData move constructor" << std::endl;
    }

    ReadData &operator=(ReadData &&other) {
//        std::cout << "ReadData move assignment" << std::endl;
        key = other.key;
        return *this;
    }

    ReadData(const ReadData &other) = default;

    ReadData &operator=(const ReadData &other) = default;
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


struct GenotypeProperties{
    std::string bases;
    uint16_t ways;
    uint16_t mat_index;
};

typedef std::vector<GenotypeProperties> GenotypeVector;

struct MutationDescription{
    uint16_t mutant_line;
    std::string from_genotype;
    std::string to_genotype;
    double line_prob;
    double genotype_prob;
    double lik;
};



typedef std::vector< std::string > SampleNames;
typedef std::unordered_map<std::string, uint32_t> SampleMap;

typedef Eigen::ArrayXd GenotypeProbs;
typedef Eigen::Array4d HaploidProbs;
typedef Eigen::Array<double, 16, 1> DiploidProbs;

typedef Eigen::Matrix4d TransitionMatrix;

typedef Eigen::ArrayXXd MutationMatrix; //array of dynamic dimensions

typedef std::array<double, 10> DiploidProbsIndex10;
typedef std::array<double, 10> Array10D;
typedef std::array<double, 4> Array4D;




#endif //ACCUMULATE_DATA_STRUCT_H
