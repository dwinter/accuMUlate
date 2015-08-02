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

namespace BoostUtils {
    using namespace std;

//    namespace po = boost::program_options;

    void ParseCommandLinkeInput(int argc, char **argv, boost::program_options::variables_map &vm);

    void ExtractInputVariables(boost::program_options::variables_map &vm, GenomeData &genome_data,
            BamTools::BamReader &experiment, BamTools::RefVector &references,
            BamTools::SamHeader &header, BamTools::Fasta &reference_genome);

    ModelParams CreateModelParams(boost::program_options::variables_map variables_map);
    
    void ValidateNfreqs(boost::any& v, const vector<string>& values,  vector<double>*)

    void CheckArgs(boost::program_options::variables_map vm)


}

#endif //BOOST_INPUT_UTILS_H_
