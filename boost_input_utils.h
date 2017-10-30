/*
 * boost_input_utils.h
 *
 *  Created on: 12/12/14
 *      Author: Steven Wu
 */


#pragma once
#ifndef BOOST_INPUT_UTILS_H_
#define BOOST_INPUT_UTILS_H_

#include <sys/stat.h>
#include <time.h>

#include <boost/program_options/variables_map.hpp>
#include <boost/program_options.hpp>
#include <api/BamReader.h>

#include "parsers.h"
#include "src/io_data/local_bamtools/bamtools_pileup_engine.h"
#include "src/io_data/local_bamtools/bamtools_fasta.h"
//#include "utils/bamtools_pileup_engine.h"
//#include "utils/bamtools_fasta.h"
#include "model.h"


struct nfreqs{
    std::vector<double> freqs;
};

void validate(boost::any& v, const vector<string>& values, nfreqs target_type, int);

namespace BoostUtils {
    using namespace std;

//    namespace po = boost::program_options;
//    static bool file_exists(const std::string &name);

    void ParseCommandLineInput(int argc, char **argv, boost::program_options::variables_map &vm);
    void ParseDenominateCommandline(int argc, char **argv, boost::program_options::variables_map &vm);
    void ExtractInputVariables(boost::program_options::variables_map &vm,
            BamTools::BamReader &experiment, BamTools::RefVector &references,
            BamTools::SamHeader &header, LocalBamToolsUtils::Fasta &reference_genome);

    ModelParams CreateModelParams(boost::program_options::variables_map variables_map);
    SampleMap ParseSamples(boost::program_options::variables_map &vm, BamTools::SamHeader &header);
//    void check_args(boost::program_options::variables_map &vm);
        
    


}

#endif //BOOST_INPUT_UTILS_H_
