#include <iostream>
#include <stdint.h>
#include <map>
#include <vector>
#include <string>

#include "boost/program_options.hpp"
#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"
#include "utils/bamtools_fasta.h"

#include "model.h"

using namespace std;
using namespace BamTools;

typedef vector< string > SampleNames;

string get_sample(string& tag){
    string res;
    for(size_t i = 0; tag[i] != '_'; i++) {
        res.push_back(tag[i]);
    }
    return(res);            
}

int find_sample_index(string s, SampleNames sv){
    for (size_t i=0; i < sv.size(); i++){
        if(s.compare(sv[i])==0){
            return i;
        }
    }
    cerr << "didn't find name " << s << endl;
    return(13); //TODO refactor this to  update sample in place
}

uint16_t base_index(char b){
    switch(b){//TODO best 'null/npos' result for a short int?
        case 'A':
        case 'a':    
            return 0;
         case 'T':
         case 't':
             return 3;
         case 'C':
         case 'c':
             return 1;
         case 'G':
         case 'g':
             return 2;
         case '-':
         case 'N':
             return 4;
         default:
             cerr << "Don't know what to make of base" << b <<endl;
             return 4;
    };
}

                             
class VariantVisitor : public PileupVisitor{
    public:
        VariantVisitor(const RefVector& bam_references, 
                       const Fasta& idx_ref,
                       ostream *out_stream,
                       SampleNames samples, 
                       const ModelParams& p,  
                       BamAlignment& ali, 
                       int q,
                       double prob_cutoff):

            PileupVisitor(), m_idx_ref(idx_ref), m_bam_ref(bam_references), 
                             m_samples(samples), m_q(q), m_params(p), m_ali(ali), m_ostream(out_stream)
            {
                nsamp = m_samples.size();
            }
        ~VariantVisitor(void) { }
    public:
         void Visit(const PileupPosition& pileupData) {
             string chr = m_bam_ref[pileupData.RefId].RefName;
             uint64_t pos  = pileupData.Position;
             m_idx_ref.GetBase(pileupData.RefId, pos, current_base);
             ReadDataVector bcalls (nsamp, ReadData{{ 0,0,0,0 }}); //fill constructor
             string tag_id;
             for(auto it = begin(pileupData.PileupAlignments);
                      it != end(pileupData.PileupAlignments); 
                      ++it){
                 int const *pos = &it->PositionInAlignment;
                 if (it->Alignment.Qualities[*pos] - 33 > m_q){
                     it->Alignment.GetTag("RG", tag_id);
                     int sindex = find_sample_index(get_sample(tag_id), m_samples);
                     uint16_t bindex  = base_index(it->Alignment.AlignedBases[*pos]);
                     if (bindex < 4){
                         bcalls[sindex].reads[bindex] += 1;
                     }
                 }
            }
            ModelInput d = {base_index(current_base), bcalls};
            double prob = TetMAProbability(m_params,d);
            if(prob >= 0.0){//NOTE: Probablity cut off is hard-coded atm
             *m_ostream << chr << '\t' << pos << '\t' << current_base << '\t' << 
                 prob << '\t' << TetMAProbOneMutation(m_params,d) << endl;          
            }
        }
    private:
        Fasta m_idx_ref; 
        RefVector m_bam_ref;
        SampleNames m_samples;
        int m_q;
        int nsamp;
        ModelParams m_params;
        BamAlignment& m_ali;
        char current_base;
        uint64_t chr_index;
        ostream* m_ostream;
            
};



int main(int argc, char** argv){

    namespace po = boost::program_options;
    string ref_file;
    po::options_description cmd("Command line options");
    cmd.add_options()
        ("help,h", "Print a help message")
        ("bam,b", po::value<string>(), "Path to BAM file")
        ("reference,r", po::value<string>(&ref_file), "Path to reference genome")
//        ("ancestor,a", po::value<string>(&anc_tag), "Ancestor RG sample ID")
        ("sample-name,s", po::value<vector <string> >(), "Sample tags")
        ("qual,q", po::value<int>()->default_value(13), 
                   "Base quality cuttoff (default = 13)")
        ("prob,p", po::value<double>()->default_value(0.1),
                   "Mutaton probability cut-off (default = 0.1)")
        ("out,o", po::value<string>()->default_value("accuMUlate_"), 
                   "Stem for output file names (default = 'accMUlate_'")
        ("config,c", po::value<string>(), "Path to config file")
        ("theta", po::value<double>(), "theta")            
        ("nfreqs", po::value<vector<double> >()->multitoken(), "")     
        ("mu", po::value<double>(), "")  
        ("seq-error", po::value<double>(), "") 
        ("phi-haploid",     po::value<double>(), "") 
        ("phi-diploid",     po::value<double>(), ""); 

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmd), vm);
    vm.notify();

    if (vm.count("help")){
        cout << cmd << endl;
        return 0;
    }

    if (vm.count("config")){
        ifstream config_stream (vm["config"].as<string>());
        po::store(po::parse_config_file(config_stream, cmd, false), vm);
        vm.notify();
    }

    ModelParams params = {
        vm["theta"].as<double>(),
        vm["nfreqs"].as<vector< double> >(),
        vm["mu"].as<double>(),
        vm["seq-error"].as<double>(), 
        vm["phi-haploid"].as<double>(), 
        vm["phi-diploid"].as<double>(),
    };

    ofstream result_stream (vm["out_name"].as<string>());

    BamReader experiment; 
    experiment.Open(vm["bam"].as<string>());
    RefVector references = experiment.GetReferenceData(); 
    Fasta reference_genome; // BamTools::Fasta
    reference_genome.Open(ref_file);
    reference_genome.CreateIndex(ref_file + ".fai");

    BamAlignment ali;
    PileupEngine pileup;
    VariantVisitor *v = new VariantVisitor(
            references,
            reference_genome, 
            &result_stream,
            vm["sample-name"].as<vector< string> >(),
            params, 
            ali, 
            vm["qual"].as<int>(), 
            vm["prob"].as<int>() );
    pileup.AddVisitor(v);
    while( experiment.GetNextAlignment(ali)){
        pileup.AddAlignment(ali);
         
    };
    pileup.Flush();
    return 0;
}
        
        



