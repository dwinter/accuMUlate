/*
 * boost_input_utils.h
 *
 *  Created on: 12/12/14
 *      Author: Steven Wu
 */


#pragma once
#ifndef BOOST_INPUT_UTILS_H_
#define BOOST_INPUT_UTILS_H_


#include <boost/program_options/variables_map.hpp>

#include <boost/program_options.hpp>
#include <time.h>


#include "parsers.h"
#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"
#include "utils/bamtools_fasta.h"

#include "model.h"
#include <sys/stat.h>

struct nfreqs{
    vector<double> freqs;
};
//    double fA,fC,fG,fT;
//    nfreqs(vector<double> pi){
//        fA = pi[0];
//        fC = pi[1];
//        fG = pi[2];
//        fT = pi[3];
//    }
//    vector<double> asVector(){
//        vector<double> v {fA, fC, fG, fT};
//        return v;
//    }
//};

namespace BoostUtils {
    using namespace std;

//    namespace po = boost::program_options;

    void ParseCommandLineInput(int argc, char **argv, boost::program_options::variables_map &vm);

    void ExtractInputVariables(boost::program_options::variables_map &vm,
            BamTools::BamReader &experiment, BamTools::RefVector &references,
            BamTools::SamHeader &header, BamTools::Fasta &reference_genome);

    ModelParams CreateModelParams(boost::program_options::variables_map variables_map);
    SampleMap ParseSamples(boost::program_options::variables_map &vm, BamTools::SamHeader &header);
    
    void ValidateNfreqs(boost::any& v, const vector<string>& values,  vector<double>*);

    void CheckArgs(boost::program_options::variables_map vm);


}

#endif //BOOST_INPUT_UTILS_H_
