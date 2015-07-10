#include <iostream>
#include <stdint.h>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <sys/stat.h>


#include "boost/program_options.hpp"
#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"
#include "utils/bamtools_fasta.h"

#include "model.h"
#include "parsers.h"

using namespace std;
using namespace BamTools;




bool include_sample(const ModelParams &params, const ReadDataVector fwd, const ReadDataVector rev,  const ReadDataVector site_data, int sindex, uint16_t ref_base, double pcut, bool central){
    
    //Can't be included if you don't have 3fwd, 3rev so check that before we do
    //any number crunching    
    //if(central == true){
        auto Fm = max_element(fwd[sindex].reads, fwd[sindex].reads+4);
        if (*Fm < 3){
            return false;
        }
        
        auto Rm = max_element(rev[sindex].reads, rev[sindex].reads+4);
        if (*Rm < 3){
            return false;
        }
   // }
    ReadDataVector _site_data = site_data;
    for(size_t i = 1; i < 4; i++){
        rotate( begin(_site_data[sindex].reads), begin(_site_data[sindex].reads) + i, end(_site_data[sindex].reads) );
        ModelInput d = { ref_base, _site_data };
        double p  = TetMAProbability(params, d); 
        if (p > 0.1){
            return true;
        }
    }
    return false;
}





void call_ancestor(const ModelParams &params, int ref_allele, const ReadData &d){
    uint16_t result[2];
    if( (d.reads[0] + d.reads[1] + d.reads[2] + d.reads[3]) == 0){
        result[0] = ref_allele;
        result[1] = ref_allele;
    }
    else{
    	DiploidProbs genotypes = DiploidSequencing(params, ref_allele, d);
        Eigen::Array33d::Index idx;
        //std::cerr << genotypes.maxCoeff() << std::endl;
        genotypes.maxCoeff(&idx);
        result[0] = idx / 4;
        result[1] = idx % 4;
    }
    std::cerr << result[0] << '\t' << result[1] << '\t';
}



    

                             
class VariantVisitor : public PileupVisitor{
    public:
        VariantVisitor(const RefVector& bam_references, 
                       const SamHeader& header,
                       const Fasta& idx_ref,
                       const SampleMap& samples, 
                       int nsamples,
                       BamAlignment& ali, 
                       int qual_cut,
                       int mapping_cut,
                       ReadDataVector &denoms,
                       ModelParams& params):

            PileupVisitor(), m_idx_ref(idx_ref), m_bam_ref(bam_references), 
                             m_header(header), m_samples(samples),m_nsamp(nsamples), 
                             m_qual_cut(qual_cut), m_ali(ali), 
                             m_denoms(denoms),
                             m_mapping_cut(mapping_cut), m_params(params)
                              { }

        ~VariantVisitor(void) { }
    public:
         void Visit(const PileupPosition& pileupData) {
             uint64_t pos  = pileupData.Position;
             uint32_t dist_to_end  = ( (pos < 500) ? pos :  (m_bam_ref[pileupData.RefId].RefLength - pos));
             bool central = dist_to_end > 500;
             m_idx_ref.GetBase(pileupData.RefId, pos, current_base);
             ReadDataVector fwd_calls (m_samples.size(), ReadData{{ 0,0,0,0 }}); 
             ReadDataVector rev_calls (m_samples.size(), ReadData{{ 0,0,0,0 }});
             for(auto it = begin(pileupData.PileupAlignments);
                      it !=  end(pileupData.PileupAlignments); 
                      ++it){
                 if( include_site(*it, m_mapping_cut, m_qual_cut) ){
                    it->Alignment.GetTag("RG", tag_id);
                    uint32_t sindex = m_samples[tag_id]; //TODO check samples existed! 
                    uint16_t bindex  = base_index(it->Alignment.QueryBases[it->PositionInAlignment]);
                    if (bindex < 4 ){
                        if(it->Alignment.IsReverseStrand() ){
                            fwd_calls[sindex].reads[bindex] += 1;
                        }
                        else{
                            rev_calls[sindex].reads[bindex]+= 1;
                        }
                    }
                }
            }
             
            uint16_t ref_base_idx = base_index(current_base);
            if (ref_base_idx < 4  ){ //TODO Model for bases at which reference is 'N' 
                ReadDataVector all_calls (m_samples.size(), ReadData{{ 0,0,0,0 }});
                for(size_t i = 0; i < m_samples.size(); i ++){
                    for(size_t j =0; j < 4; j++){
                        uint16_t sum_calls = fwd_calls[i].reads[j] + rev_calls[i].reads[j];
                        all_calls[i].reads[j] = sum_calls;
                    }
                }

                for(size_t i  = 1; i < m_samples.size(); i++){
                    if( include_sample(m_params, fwd_calls, rev_calls, all_calls, i, ref_base_idx, 0.1, central) ){
                        m_denoms[i].reads[ref_base_idx] += 1;
                    }
                }
            }
         }

         






                 
//                ModelInput d = {ref_base_idx, bcalls};
//                double prob_one = TetMAProbOneMutation(m_params,d);
//                double prob = TetMAProbability(m_params, d);
//                if(prob >= m_prob_cut){
//                     *m_ostream << m_bam_ref[pileupData.RefId].RefName << '\t'
//                                << pos << '\t' 
//                                << current_base << '\t' 
//                                << prob << '\t' 
//                                << prob_one << '\t' 
//                                << endl;          
//                }
//            }
    private:
        RefVector m_bam_ref;
        SamHeader m_header;
        Fasta m_idx_ref; 
        SampleMap m_samples;
        int m_nsamp;
        BamAlignment& m_ali;
        int m_qual_cut;
        int m_mapping_cut;
        char current_base;
        string tag_id;
        uint64_t chr_index;
        ReadDataVector& m_denoms;
        ModelParams m_params;
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
     
        ("intervals,i", po::value<string>(), "Path to bed file");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmd), vm);

    if (vm.count("help")){
        cout << cmd << endl;
        return 0;
    }
    ModelParams params = {
        0.0001,
        {0.388, 0.112, 0.112, 0.338},
        1e-8,
        0.01,
        0.01, 
        0.005
    };

    vm.notify();
    string bam_path = vm["bam"].as<string>();
    string index_path = vm["bam-index"].as<string>();
    if(index_path == ""){
        index_path = bam_path + ".bai";
    }   


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
        reference_genome.CreateIndex(faidx_path);
    }
    reference_genome.Open(ref_file, faidx_path);

    // Map readgroups to samples
    // TODO: this presumes first sample is ancestor. True for our data, not for
    // others.
    // First map all sample names to an index for ReadDataVectors
    SampleMap name_map;
    uint16_t sindex = 0;
    for(auto it = header.ReadGroups.Begin(); it!= header.ReadGroups.End(); it++){
        if(it->HasSample()){
            auto s  = name_map.find(it->Sample);
            if( s == name_map.end()){ // not in there yet
                name_map[it->Sample] = sindex;
                sindex += 1;
            }
        }
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

    ReadDataVector denoms (sindex, {0,0,0,0} );

    VariantVisitor *v = new VariantVisitor(
            references,
            header,
            reference_genome, 
//            vm["sample-name"].as<vector< string> >(),
            samples,
            sindex,
            ali, 
            vm["qual"].as<int>(), 
            vm["mapping-qual"].as<int>(),
            denoms,
            params            
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
    for(size_t i = 0; i < sindex; i++){
        for( size_t j = 0; j < 4; j++){
            cout << denoms[i].reads[j] << '\t'; 
        }
    }
    cout << endl;
    return 0;
}


