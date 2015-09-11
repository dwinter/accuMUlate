#include <iostream>
#include <stdint.h>
#include <map>
#include <vector>
#include <string>
#include <sys/stat.h>
#include <chrono>


#include "boost/program_options.hpp"
#include "api/BamReader.h"


#include "src/io_data/local_bamtools/bamtools_pileup_engine.h"
#include "src/io_data/local_bamtools/bamtools_fasta.h"
//#include "utils/bamtools_pileup_engine.h"
//#include "utils/bamtools_fasta.h"


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
            params,
            ali,
            vm["qual"].as<int>(),
            vm["mapping-qual"].as<int>(),
            vm["prob"].as<double>()
        );

    pileup.AddVisitor(v);

    std::cerr.setstate(std::ios_base::failbit) ; //Supress bamtool camplain, "Pileup::Run() : Data not sorted correctly!"

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    if (vm.count("intervals")){
        BedFile bed (vm["intervals"].as<string>());
        BedInterval region;
        while(bed.get_interval(region) == 0){
            int ref_id = experiment.GetReferenceID(region.chr);
            experiment.SetRegion(ref_id, region.start, ref_id, region.end);
            while( experiment.GetNextAlignment(ali) ){
                pileup.AddAlignment(ali);
            }
        }
    }
    else{
        while( experiment.GetNextAlignment(ali)){
            pileup.AddAlignment(ali);
        }
    }
    pileup.Flush();

    std::cerr.clear() ; // Add std::cerr back

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t start_time = std::chrono::system_clock::to_time_t(start);
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "Started at "<< std::ctime(&start_time) << "EM finished at " << std::ctime(&end_time)  << "\nElapsed time: " << elapsed_seconds.count() << "s\n";

    return 0;
}


