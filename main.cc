#include <iostream>
#include <stdint.h>
#include <map>
#include <vector>
#include <string>

#include <boost/program_options.hpp>
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
                       GenomeData& all_the_data,
//                      const ModelParams& p,  
                       const SampleMap& samples, 
                       BamAlignment& ali, 
                       int qual_cut,
                       int mapping_cut,
                       double prob_cut):

            PileupVisitor(), m_idx_ref(idx_ref), m_bam_ref(bam_references), 
                             m_header(header), m_samples(samples), 
                             m_qual_cut(qual_cut), m_ali(ali), 
                             m_all_the_data(all_the_data), m_prob_cut(prob_cut),
                             m_mapping_cut(mapping_cut)
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
                      it !=  end(pileupData.PileupAlignments); 
                      ++it){
                 if( include_site(*it, m_mapping_cut, m_qual_cut) ){
                    it->Alignment.GetTag("RG", tag_id);
                    string sm =  m_header.ReadGroups[tag_id].Sample;
                    uint32_t sindex = m_samples[sm]; //TODO check samples existed! 
                    uint16_t bindex  = base_index(it->Alignment.QueryBases[it->PositionInAlignment]);
                    if (bindex < 4 ){
                        bcalls[sindex].reads[bindex] += 1;
                    }
                }
            }
            uint16_t ref_base_idx = base_index(current_base);
            if (ref_base_idx < 4  ){ //TODO Model for bases at which reference is 'N'
//                ModelInput d = {ref_base_idx, bcalls};
                m_all_the_data.push_back(ModelInput{ ref_base_idx, bcalls });
//                ModelInput d = {ref_base_idx, bcalls};
//                double prob_one = TetMAProbOneMutation(m_params,d);
//                double prob = TetMAProbability(m_params, d);
//                if(prob >= m_prob_cut){
//                     *m_ostream << chr << '\t' 
//                                << pos << '\t' 
//                                << current_base << '\t' 
//                                << prob << '\t' 
//                                << prob_one << '\t' 
//                               << endl;          
//                }
            }
         }
    private:
        RefVector m_bam_ref;
        SamHeader m_header;
        Fasta m_idx_ref; 
        GenomeData& m_all_the_data;
        SampleMap m_samples;
        BamAlignment& m_ali;
//        ModelParams m_params;
        int m_qual_cut;
        int m_mapping_cut;
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
        ("bam-index,x", po::value<string>()->default_value(""), "Path to BAM index, (defalult is <bam_path>.bai")
        ("reference,r", po::value<string>(&ref_file)->required(),  "Path to reference genome")
//       ("ancestor,a", po::value<string>(&anc_tag), "Ancestor RG sample ID")
//        ("sample-name,s", po::value<vector <string> >()->required(), "Sample tags")
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
//    ModelParams params = {
//        vm["theta"].as<double>(),
//        vm["nfreqs"].as<vector< double> >(),
//        vm["mu"].as<double>(),
//        vm["seq-error"].as<double>(), 
//        vm["phi-haploid"].as<double>(), 
//        vm["phi-diploid"].as<double>(),
//    };
    string bam_path = vm["bam"].as<string>();
    string index_path = vm["bam-index"].as<string>();
    if(index_path == ""){
        index_path = bam_path + ".bai";
    }   

//    ofstream result_stream (vm["out"].as<string>());
    //TODO: check sucsess of all these opens/reads:
    BamReader experiment; 
    experiment.Open(bam_path);
    experiment.OpenIndex(index_path);
    RefVector references = experiment.GetReferenceData(); 
    SamHeader header = experiment.GetHeader();
    Fasta reference_genome; // BamTools::Fasef_file);
    reference_genome.Open(ref_file, ref_file+ ".fai");
//    reference_genome.CreateIndex(ref_file + ".fai");
    PileupEngine pileup;
    BamAlignment ali;


//  Assign some memory for the big list
    uint32_t total_len = 0;
    if (vm.count("intervals")){
        BedFile bed (vm["intervals"].as<string>());
        BedInterval region;
        while(bed.get_interval(region) == 0){
            total_len += (region.end - region.start);
        }
    }
    else {
        
        for_each(references.begin(),references.end(),[&](RefData chrom){
            total_len += chrom.RefLength;
         });
    }
    GenomeData base_counts;
    base_counts.reserve(total_len);
    
    SampleMap samples;
    uint16_t sindex = 0;
    for(auto it = header.ReadGroups.Begin(); it!= header.ReadGroups.End(); it++){
        if(it->HasSample()){
            auto s  = samples.find(it->Sample);
            if( s == samples.end()){ // not in there yet
                samples[it->Sample] = sindex;
                sindex += 1;
            }
        }
    }

    VariantVisitor *v = new VariantVisitor(
            references,
            header,
            reference_genome, 
            base_counts,
//            &result_stream,
            samples,
//            params, 
            ali, 
            vm["qual"].as<int>(), 
            vm["mapping-qual"].as<int>(),
            vm["prob"].as<double>()
        );
    pileup.AddVisitor(v);
//  TODO: Only allocate interval-sized memory vector   
//  if intervals are set
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
        clock_t t;
        uint64_t ali_counter = 0;
        t = clock();
        BamAlignment ali;
        while( experiment.GetNextAlignment(ali)){
            pileup.AddAlignment(ali);
            ali_counter += 1;
            if (ali_counter % 1000000 == 0){
                t = clock() - t;
                cout << "Processed 1 million reads (" 
                     << ((float)t)/CLOCKS_PER_SEC
                     << " seconds)" << endl;
            }
        }
    } 
    pileup.Flush();
    cout << base_counts.size() << endl;
    return 0;
}


