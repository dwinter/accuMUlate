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
#include "parsers.h"

using namespace std;
using namespace BamTools;

                             
class VariantVisitor : public PileupVisitor{
    public:
        VariantVisitor(const RefVector& bam_references, 
                       const SamHeader& header,
                       const Fasta& idx_ref,
                       ostream *out_stream,
                       SampleNames samples, 
                       const ModelParams& p,  
                       BamAlignment& ali, 
                       int qual_cut,
                       double prob_cut):

            PileupVisitor(), m_idx_ref(idx_ref), m_bam_ref(bam_references), 
                             m_header(header), m_samples(samples), 
                             m_qual_cut(qual_cut), m_params(p), m_ali(ali), 
                             m_ostream(out_stream), m_prob_cut(prob_cut)
                              { }
        ~VariantVisitor(void) { }
    public:
         void Visit(const PileupPosition& pileupData) {
             string chr = m_bam_ref[pileupData.RefId].RefName;
             uint64_t pos  = pileupData.Position;
             m_idx_ref.GetBase(pileupData.RefId, pos, current_base);
             ReadDataVector bcalls (m_samples.size(), ReadData{{ 0,0,0,0 }}); 
             string tag_id;
             for(auto it = begin(pileupData.PileupAlignments);
                      it != end(pileupData.PileupAlignments); 
                      ++it){
                 int const *pos = &it->PositionInAlignment;
                 if (it->Alignment.Qualities[*pos] - 33 >= m_qual_cut){
                     it->Alignment.GetTag("RG", tag_id);
                     string sm =  m_header.ReadGroups[tag_id].Sample;
                     uint32_t sindex = find_sample_index(sm, m_samples); //TODO check samples existed! 
                     uint16_t bindex  = base_index(it->Alignment.AlignedBases[*pos]);
                     if (bindex < 4 ){
                         bcalls[sindex].reads[bindex] += 1;
                     }
                 }
            }
            uint16_t ref_base_idx = base_index(current_base);
            if (ref_base_idx < 4  ){ //TODO Model for bases at which reference is 'N'
                ModelInput d = {ref_base_idx, bcalls};
                double prob_one = TetMAProbOneMutation(m_params,d);
                double prob = TetMAProbability(m_params, d);
                if(prob >= m_prob_cut){
                     *m_ostream << chr << '\t' 
                                << pos << '\t' 
                                << current_base << '\t' 
                                << prob << '\t' 
                                << prob_one << '\t' 
                                << endl;          
                }
            }
        }
    private:
        RefVector m_bam_ref;
        SamHeader m_header;
        Fasta m_idx_ref; 
        ostream* m_ostream;
        SampleNames m_samples;
        BamAlignment& m_ali;
        ModelParams m_params;
        int m_qual_cut;
        double m_prob_cut;
        char current_base;
        uint64_t chr_index;
            
};



int main(int argc, char** argv){

    namespace po = boost::program_options;
    string ref_file;
    string config_path;
    po::options_description cmd("Command line options");
    cmd.add_options()
        ("help,h", "Print a help message")
        ("bam,b", po::value<string>()->required(), "Path to BAM file")
        ("reference,r", po::value<string>(&ref_file)->required(),
                        "Path to reference genome")
//        ("ancestor,a", po::value<string>(&anc_tag), "Ancestor RG sample ID")
        ("sample-name,s", po::value<vector <string> >()->required(), "Sample tags")
        ("qual,q", po::value<int>()->default_value(13), 
                   "Base quality cuttoff")
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
    };

    ofstream result_stream (vm["out"].as<string>());

    BamReader experiment; 
    experiment.Open(vm["bam"].as<string>());
    RefVector references = experiment.GetReferenceData(); 
    SamHeader header = experiment.GetHeader();
    Fasta reference_genome; // BamTools::Fasta
    reference_genome.Open(ref_file);
    reference_genome.CreateIndex(ref_file + ".fai");
    PileupEngine pileup;
    BamAlignment ali;
    VariantVisitor *v = new VariantVisitor(
            references,
            header,
            reference_genome, 
            &result_stream,
            vm["sample-name"].as<vector< string> >(),
            params, 
            ali, 
            vm["qual"].as<int>(), 
            vm["prob"].as<double>() );
    pileup.AddVisitor(v);
   
    if (vm.count("intervals")){
        BedFile bed (vm["intervals"].as<string>());
        FastaReference my_ref (ref_file + ".fai");
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

        BamAlignment ali;
        while( experiment.GetNextAlignment(ali)){
            pileup.AddAlignment(ali);
        };  
    pileup.Flush();
    return 0;
    }
}
        
        



