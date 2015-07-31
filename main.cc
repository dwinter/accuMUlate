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
                       MutationMatrix prob_paths_m,
                       MutationMatrix prob_paths_nm):

        ReadDataVisitor(bam_references, idx_ref, samples, p, ali, qual_cut, mapping_cut),
                        m_ostream(out_stream), m_prob_cut(prob_cut),
                        m_prob_paths_m(prob_paths_m), m_prob_paths_nm(prob_paths_nm)
                        {}
        
        ~VariantVisitor(void) { }
    public:
         void Visit(const PileupPosition& pileupData) {
            if (GatherReadData(pileupData) ){
                double prob = TetMAProbability(m_params, site_data, m_prob_paths_m, m_prob_paths_nm);
                if(prob >= m_prob_cut){
//                    double prob_one = TetMAProbOneMutation(m_params, site_data);
                    double prob_one = 0.0;
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
        MutationMatrix m_prob_paths_m;
        MutationMatrix m_prob_paths_nm;
        ostream* m_ostream;
        double m_prob_cut;

};



int main(int argc, char** argv){

    namespace po = boost::program_options;
    string ref_file;
    string config_path;
    string anc_tag;
    po::options_description cmd("Command line options");
    cmd.add_options()
        ("help,h", "Print a help message")
        ("bam,b", po::value<string>()->required(), "Path to BAM file")
        ("bam-index,x", po::value<string>()->default_value(""), "Path to BAM index, (defalult is <bam_path>.bai")
        ("reference,r", po::value<string>(&ref_file)->required(),  "Path to reference genome")
        ("ancestor,a", po::value<string>(&anc_tag), "Ancestor RG sample ID")
        ("sample-name,s", po::value<vector <string> >()->required(), "Sample tags to include")
        ("qual,q", po::value<int>()->default_value(13), 
                   "Base quality cuttoff")
        
        ("mapping-qual,m", po::value<int>()->default_value(13), 
                    "Mapping quality cuttoff")
     
        ("prob,p", po::value<double>()->default_value(0.1),
                   "Mutaton probability cut-off")
        ("out,o", po::value<string>()->default_value("acuMUlate_result.tsv"),
                    "Out file name")
        ("intervals,i", po::value<string>(), "Path to bed file")
        ("config,c", po::value<string>(), "Path to config file")
        ("theta", po::value<double>()->required(), "theta")            
        ("nfreqs", po::value<vector<double> >()->multitoken(), "")     
        ("mu", po::value<double>()->required(), "")  
        ("seq-error", po::value<double>()->required(), "")
        ("ploidy-ancestor", po::value<int>()->default_value(2), "")
        ("ploidy-descendant", po::value<int>()->default_value(2), "")
        ("phi-haploid",     po::value<double>()->required(), "") 
        ("phi-diploid",     po::value<double>()->required(), "");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmd), vm);

    if (vm.count("help")){
        cout << cmd << endl;
        return 0;
    }

    if (vm.count("config")){
        ifstream config_stream (vm["config"].as<string>());
        po::store(po::parse_config_file(config_stream, cmd, false), vm);
    }

    vm.notify();
    ModelParams params = {
        vm["theta"].as<double>(),
        vm["nfreqs"].as<vector< double> >(),
        vm["mu"].as<double>(),
        vm["seq-error"].as<double>(), 
        vm["phi-haploid"].as<double>(), 
        vm["phi-diploid"].as<double>(),
        vm["ploidy-ancestor"].as<int>(),
        vm["ploidy-descendant"].as<int>()
    };
    string bam_path = vm["bam"].as<string>();
    string index_path = vm["bam-index"].as<string>();
    if(index_path == ""){
        index_path = bam_path + ".bai";
    }   

    ofstream result_stream (vm["out"].as<string>());
    // Start setiing up files
    //TODO: check sucsess of all these opens/reads:

    BamReader experiment; 
    experiment.Open(bam_path);
    experiment.OpenIndex(index_path);
    RefVector references = experiment.GetReferenceData(); 
    SamHeader header = experiment.GetHeader();

    
    //Fasta reference
    Fasta reference_genome; // BamTools::Fasta
    struct stat file_info;
    string faidx_path = ref_file + ".fai";
    if (stat(faidx_path.c_str(), &file_info) != 0){
        reference_genome.Open(ref_file);
        reference_genome.CreateIndex(faidx_path);
    }
    reference_genome.Open(ref_file, faidx_path);

    // Map readgroups to samples
    // First check it we want to use a given sample 
    // map all 'keeper' sample names to an index for ReadDataVectors use
    // (unsigned) -1 as a "missing" code
    // 
    SampleMap name_map;
    vector<string> keepers =  vm["sample-name"].as< vector<string> >();
    uint16_t sindex = 1;
    bool ancestor_in_BAM = false;
    for(auto it = header.ReadGroups.Begin(); it!= header.ReadGroups.End(); it++){
        if(it->HasSample()){
            if (find(keepers.begin(), keepers.end(), it->Sample) == keepers.end()){
                if(it->Sample == anc_tag ){    
                      name_map[it->Sample] = 0; 
                      ancestor_in_BAM = true;
                } else {      
                    name_map[it->Sample] = std::numeric_limits<uint32_t>::max()  ;           
                    cerr << "Warning: excluding data from '" << it->Sample <<
                         "' which is included in the BAM file but not the list of included samples" << endl;
                }
            }else {
                auto s  = name_map.find(it->Sample);
                if( s == name_map.end()){ 
                    name_map[it->Sample] = sindex;
                    sindex += 1;
                }
            }
        }
    }
    if(!ancestor_in_BAM){
        cerr << "Error: No data for ancestral sample '" << anc_tag << 
            "' in the specifified BAM file. Check the sample tags match" << endl;
        exit(1);
    }
    // And now, go back over the read groups to map RG:sample index
    SampleMap samples;
    for(auto it = header.ReadGroups.Begin(); it!= header.ReadGroups.End(); it++){
        if(it->HasSample()){
            samples[it->ID] = name_map[it->Sample];  
        }
    }

    PileupEngine pileup;
    BamAlignment ali;

    MutationMatrix mt = MutationAccumulation(params, true);
    MutationMatrix m = MutationAccumulation(params, false);
    MutationMatrix nm = m - mt;
    cerr << mt <<endl;

    VariantVisitor *v = new VariantVisitor(
            references,
            reference_genome, 
            &result_stream,
//            vm["sample-name"].as<vector< string> >(),
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


