#include <iostream>
#include <stdint.h>
#include <string>
#include <algorithm>
#include <vector>
#include <numeric> 

#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"
#include "boost/program_options.hpp"

#include "model.h"
#include "parsers.h"

using namespace std;
using namespace BamTools;


class SampleSiteData{

    public:
        vector<uint16_t> BQs;
        vector<uint16_t> MQs;
        uint32_t fwd_reads;
        uint32_t rev_reads;
        ReadData base_calls;

        SampleSiteData(){
          base_calls.key = 0;
          fwd_reads = 0;
          rev_reads = 0;
        }
                  
        int get_genotype(){
           //TODO: call genotypes from the model
           //These have already been called for mutation-ness, and are haploid
           //so, to make a start, we are just calling the most common base         
           if (base_calls.key == 0){
              return -1;
           } 
           return(distance(base_calls.reads, max_element(base_calls.reads, base_calls.reads + 4)));
        }

        void import_alignment(const BamAlignment& al, const int& pos, const int& bindex){
            uint16_t b_index = base_index(al.AlignedBases[pos]);
            base_calls.reads[bindex] += 1;
            BQs.push_back(al.Qualities[pos]);
            MQs.push_back(al.MapQuality);
            if(al.IsReverseStrand()){ 
                rev_reads += 1; 
            } 
            else{ 
                fwd_reads += 1; 
            }
        }
};

 
class ExperimentSiteData{
    public:
        string m_initial_data;
        vector<SampleSiteData> sample_data;
        SampleNames snames;

        ExperimentSiteData(SampleNames sn, string initial_data ){
            snames = sn;
            m_initial_data = initial_data;
            //fill constructor doesn't work here?
            for(size_t i= 0; i < sn.size(); i++){
                sample_data.push_back(SampleSiteData());
            }
        }

        void summarize(ostream* out_stream){
            // call genotypes, keep track of each sample and number of each 
            // possible allele
            *out_stream << m_initial_data;
            vector<uint16_t> genotypes;
            ReadData gfreqs;
            gfreqs.key = 0;
            for(size_t i=0; i < sample_data.size(); i++){
                uint16_t g = sample_data[i].get_genotype();
                if(g < 4){//no data == npos
                    gfreqs.reads[g] += 1;
                }
                genotypes.push_back(g);             
            }
            //Now find the mutant (should be the only one with it's allele)
            uint16_t mutant_base;
            uint16_t n_mutant = 0;
            for(size_t i = 0; i<4; i++){
                if (gfreqs.reads[i] == 1){
                    n_mutant += 1;
                    mutant_base = i;
                }
            }
            if (n_mutant != 1){
                *out_stream  << "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" << endl;
                cerr << "skipping " << m_initial_data << ". Read Matrix:" << endl;
                for (size_t i=0; i < sample_data.size(); i++){
                    cerr << snames[i] << '\t';
                    for (size_t j=0; j<4; j++){
                        cerr << sample_data.base_call.reads[j] << '\t';
                    }
                }
                cerr << endl;
                return; // warn?
            }
            auto it = find_if(genotypes.begin(), genotypes.end(), 
                    [&](int v) {return v==mutant_base;});
            uint32_t mutant = distance(genotypes.begin(), it);
            //summarise the data

            int mutant_alleles = 0;
            int mutant_alleles_denom = 0;
            int wt_MQs = 0;
            int wt_MQs_denom = 0;
            auto s_mutant_MQs = accumulate (sample_data[mutant].MQs.begin(),
                                            sample_data[mutant].MQs.end(), 0);
            double mutant_MQs = s_mutant_MQs/(double)sample_data[mutant].MQs.size();
            for (size_t i=0; i < sample_data.size(); i++){
                if (i != mutant){
                    mutant_alleles += sample_data[i].base_calls.reads[mutant_base];
                    mutant_alleles_denom += sample_data[i].fwd_reads + sample_data[i].rev_reads;
                    wt_MQs += accumulate(sample_data[i].MQs.begin(), sample_data[i].MQs.end(), 0); 
                    wt_MQs_denom += sample_data[i].MQs.size();     
                }
            }
            double xbar_MQs = (double)wt_MQs/wt_MQs_denom;
            double mutant_freq = (double)mutant_alleles/mutant_alleles_denom;
            double mutant_strain_freq = double(sample_data[mutant].base_calls.reads[mutant_base]) / (sample_data[mutant].fwd_reads + sample_data[mutant].rev_reads);
            *out_stream  << mutant_base << '\t'
                          << snames[mutant] << '\t'
                          << mutant_strain_freq<< '\t'
                          << mutant_freq << '\t' 
                          << sample_data[mutant].fwd_reads << '\t'
                          << sample_data[mutant].rev_reads << '\t'
                          << mutant_MQs << '\t'
                          << xbar_MQs  << endl;
            }
            
        
};


class FilterVisitor: public PileupVisitor{
    public: 
        FilterVisitor(BamAlignment& ali, 
                      const SampleNames& samples, 
                      ostream* out_stream,
                      int ref_pos,
                      string input_data):
            PileupVisitor(), m_samples(samples), m_ref_pos(ref_pos), 
                             m_out_stream(out_stream), m_initial_data(initial_data)
            {  } 
        ~FilterVisitor(void) { }


    public:
        void Visit(const PileupPosition& pileupData){
            if (pileupData.Position != m_ref_pos){
                return;
            }

            auto target_site = ExperimentSiteData(m_samples);
            for (auto it =  pileupData.PileupAlignments.begin();
                      it != pileupData.PileupAlignments.end();
                      it++){
                int const *pos = &it->PositionInAlignment;
                if(it->Alignment.Qualities[*pos] > 46){//TODO user-defined qual cut 
                    uint16_t b_index = base_index(it->Alignment.AlignedBases[*pos]);
                    if (b_index < 4){
                        string tag_id;
                        it->Alignment.GetTag("RG", tag_id);
                        uint32_t sindex = find_sample_index(get_sample(tag_id),m_samples);
                        target_site.sample_data[sindex].import_alignment(it->Alignment, *pos, b_index);               
                    }
                }
            }
        
        target_site.summarize(m_out_stream);
    }

    private:
        SampleNames m_samples;
        ostream* m_out_stream;
        int m_ref_pos;
       // ExperimentSiteData target_site;
    
};


void test(){

    //setup
    auto all_t = SampleSiteData();
    all_t.base_calls.key = 0;
    all_t.base_calls.reads[3] =  10 ;
    all_t.fwd_reads = 5;
    all_t.rev_reads = 5;

    auto mix = SampleSiteData();
    mix.base_calls.key = 0;
    mix.base_calls.reads[1] = 5;
    mix.base_calls.reads[2] = 5;
    mix.fwd_reads = 5;
    mix.rev_reads = 5;

    auto empty_site = SampleSiteData();
    empty_site.base_calls.key = 0;

    //test
    cerr << "calling genotypes " <<  endl;
    cerr << (all_t.get_genotype() == 3) << endl;
    cerr << (mix.get_genotype() == 1) << endl;
    cerr << (empty_site.get_genotype() == -1) << endl;

    //setup
//    auto out = &cerr;
    SampleNames sn = {"test1", "test2", "test3", "test4"};
    ExperimentSiteData test_site = ExperimentSiteData(sn);
    test_site.sample_data[0].base_calls.reads[1] = 10;
    test_site.sample_data[0].fwd_reads = 2;
    test_site.sample_data[0].rev_reads = 8;
    test_site.sample_data[1].base_calls.reads[1] = 10;
    test_site.sample_data[1].fwd_reads = 5;
    test_site.sample_data[1].rev_reads = 5;
    test_site.sample_data[2].base_calls.reads[0] = 10;
    test_site.sample_data[2].fwd_reads = 10;
    test_site.sample_data[2].rev_reads = 0;

    cerr << "summarising sites... " << endl;
    test_site.summarize(&cerr);
}


int main(int argc, char* argv[]){
    string bam_path;
    string input_path;
    vector< string > sample_names;
    namespace po = boost::program_options;
    po::options_description cmd("Command line args");
    cmd.add_options()
        ("help,h", "Print a help message")
        ("bam,b", po::value<string>(&bam_path)->required(), "Path to BAM file")
        ("input,i", po::value<string>(&input_path)->required(), "Path to results file")
        ("sample-name,s", po::value<vector <string> >(&sample_names)->required(), "Sample tags")
        ("config,c", po::value<string>(), "Path to config file")
        ("out,o", po::value<string>()->default_value("filtered_result.tsv"),
                    "Out file name");

 
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmd), vm);

    if (vm.count("help")){
        cout << cmd <<endl;
        return 0;
    }

    if (vm.count("config")){
        ifstream config_stream(vm["config"].as<string>());
        po::store(po::parse_config_file(config_stream, cmd, false), vm);
    }

    vm.notify();

    ofstream outfile (vm["out"].as<string>());

    ifstream putations(input_path);
    BamReader experiment;
    experiment.Open(bam_path);
    string L;
    while(getline(putations, L)){    
    
        PileupEngine pileup;
        BamAlignment ali;

        size_t i = 0;
        string chr;
        for(; L[i] != '\t'; ++i){
            chr.push_back(L[i]);
        }
        string pos_s;
        i += 1;
        for(; L[i] !='\t'; i++){
            pos_s.push_back(L[i]);
        }
        int pos = stoul(pos_s);
        int ref_id = experiment.GetReferenceID(chr);
        cout <<  chr  << '\t'<< pos << '\t' ;
        experiment.SetRegion(ref_id, pos-1, ref_id, pos);
        FilterVisitor *f = new FilterVisitor(ali, 
                                             sample_names,
                                             &outfile,
                                             pos, L);
        pileup.AddVisitor(f);
        while( experiment.GetNextAlignment(ali) ) {
            pileup.AddAlignment(ali);
        }
        pileup.Flush();
    }
    return 0;
}




            






