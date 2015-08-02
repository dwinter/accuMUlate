#include <iostream>
#include <stdint.h>
#include <map>
#include <vector>
#include <string>
#include <sys/stat.h>



#include "boost/program_options.hpp"
#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"
#include "utils/bamtools_fasta.h"

#include "boost_input_utils.h"
#include "model.h"
#include "parsers.h"

using namespace std;
using namespace BamTools;

class VariantVisitor : public ReadDataVisitor{
    public:
        VariantVisitor(const RefVector& bam_references,
                       Fasta& idx_ref,
                       ostream *out_stream,
                       SampleMap& samples,
                       const ModelParams& p,
                       BamAlignment& ali,
                       int qual_cut,
                       int mapping_cut,
                       double prob_cut,
                       MutationMatrix mut_paths,
                       MutationMatrix nomut_paths):

        ReadDataVisitor(bam_references, idx_ref, samples, p, ali, qual_cut, mapping_cut),
                        m_ostream(out_stream), m_prob_cut(prob_cut),
                        m_mut_paths(mut_paths), m_nomut_paths(nomut_paths)
                        {}

        ~VariantVisitor(void) { }
    public:
         void Visit(const PileupPosition& pileupData) {
            if (GatherReadData(pileupData) ){
                double prob = TetMAProbability(m_params, site_data, m_mut_paths, m_nomut_paths);
                if(prob >= m_prob_cut){
                    double prob_one = TetMAProbOneMutation(m_params, site_data, m_mut_paths, m_nomut_paths);
                         *m_ostream << m_bam_references[pileupData.RefId].RefName << '\t'
                                    << pileupData.Position << '\t'
                                    << current_base << '\t'
                                    << prob << '\t'
                                    << prob_one << '\t'
                                    << endl;
                }
            }
        }
    private:
        MutationMatrix m_mut_paths;
        MutationMatrix m_nomut_paths;
        ostream* m_ostream;
        double m_prob_cut;

};



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
    Fasta reference_genome; // BamTools::Fasta

    BoostUtils::ExtractInputVariables(vm,experiment, references,header,reference_genome ); 
    ModelParams params = BoostUtils::CreateModelParams(vm);
    SampleMap samples = BoostUtils::ParseSamples(vm, header);

    PileupEngine pileup;
    BamAlignment ali;

    MutationMatrix mt = MutationAccumulation(params, true);
    MutationMatrix m = MutationAccumulation(params, false);
    MutationMatrix nm = m - mt;

    VariantVisitor *v = new VariantVisitor(
            references,
            reference_genome,
            &result_stream,
            samples,
            params,
            ali,
            vm["qual"].as<int>(),
            vm["mapping-qual"].as<int>(),
            vm["prob"].as<double>(),
            m, nm
        );
    pileup.AddVisitor(v);

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
    return 0;
}


