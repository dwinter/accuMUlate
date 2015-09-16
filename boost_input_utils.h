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
#include "src/io_data/local_bamtools/bamtools_pileup_engine.h"
#include "src/io_data/local_bamtools/bamtools_fasta.h"
//#include "utils/bamtools_pileup_engine.h"
//#include "utils/bamtools_fasta.h"

#include "model.h"
#include <sys/stat.h>

struct nfreqs{
    vector<double> freqs;
};

void validate(boost::any& v, const vector<string>& values, nfreqs target_type, int);

namespace BoostUtils {
    using namespace std;

//    namespace po = boost::program_options;

    void ParseCommandLineInput(int argc, char **argv, boost::program_options::variables_map &vm);

    void ExtractInputVariables(boost::program_options::variables_map &vm,
            BamTools::BamReader &experiment, BamTools::RefVector &references,
            BamTools::SamHeader &header, LocalBamToolsUtils::Fasta &reference_genome);

    ModelParams CreateModelParams(boost::program_options::variables_map variables_map);
    SampleMap ParseSamples(boost::program_options::variables_map &vm, BamTools::SamHeader &header);


}

#endif //BOOST_INPUT_UTILS_H_
