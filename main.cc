#include <iostream>
#include <stdint.h>
#include <map>
#include <vector>
#include <string>
#include <sys/stat.h>


#include "boost/program_options.hpp"
#include "api/BamReader.h"

#include "src/io_data/local_bamtools/bamtools_pileup_engine.h"
#include "src/io_data/local_bamtools/bamtools_fasta.h"

#include "boost_input_utils.h"
#include "variant_visitor.h"

using namespace std;
using namespace BamTools;


int main(int argc, char** argv){


    boost::program_options::variables_map vm;
    BoostUtils::ParseCommandLineInput(argc, argv, vm);
    
    
    streambuf * buf;
    ofstream of;
    string out_name =  vm["out"].as<string>();

    if(out_name == "") {
        buf = cout.rdbuf();
    } 
    else {
        of.open(out_name);
        buf = of.rdbuf();
    }

    ostream result_stream(buf);

    BamReader experiment;
    RefVector references;
    SamHeader header;
    LocalBamToolsUtils::Fasta reference_genome; // BamTools::Fasta

    BoostUtils::ExtractInputVariables(vm, experiment, references, header, reference_genome);
    ModelParams params = BoostUtils::CreateModelParams(vm);
    SampleMap samples = BoostUtils::ParseSamples(vm, header);

    LocalBamToolsUtils::PileupEngine pileup;
    BamAlignment ali;

    VariantVisitor *v = new VariantVisitor(
            references,
            reference_genome,
            &result_stream,
            samples,
            vm["sample-name"].as< vector<string> >(),
            params,
            ali,
            vm["qual"].as<int>(),
            vm["mapping-qual"].as<int>(),
            vm["prob"].as<double>()
        );

    pileup.AddVisitor(v);
                                                  


    if (vm.count("intervals")){
        BedFile bed (vm["intervals"].as<string>());
        BedInterval region;
        while(bed.get_interval(region) == 0){
            int ref_id = experiment.GetReferenceID(region.chr);
            experiment.SetRegion(ref_id, region.start, ref_id, region.end);
            v->SetRegion(region);
            while( experiment.GetNextAlignment(ali) ){   
                pileup.AddAlignment(ali);                
            }
        pileup.Flush();
        }
    }
    else{
        while( experiment.GetNextAlignment(ali)){
            pileup.AddAlignment(ali);
        }
    }
    pileup.Flush();
    delete v;

    return 0;
}


