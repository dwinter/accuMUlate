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
        //These are 16bit unsigned ints - can handle 2^16/41 ~ 1600 quality
        //scores without overflowing
        ReadData fwd_reads;
        ReadData rev_reads;
        ReadData all_reads;
        ReadData MQ;
        ReadData BQ;
        uint16_t depth;


        SampleSiteData(){
          MQ.key = 0;
          BQ.key = 0;
          fwd_reads.key = 0;
          rev_reads.key = 0;
          all_reads.key= 0;
          depth = 0;
        }
                  
        uint16_t get_genotype(){
           //TODO: call genotypes from the model
           //These have already been called for mutation-ness, and are haploid
           //so, to make a start, we are just calling the most common base 
           if (fwd_reads.key == 0 && rev_reads.key == 0){
               return -1;
           }
           for(size_t i=0; i<4; i++){
               int n = fwd_reads.reads[i] + rev_reads.reads[i];
               all_reads.reads[i] += n;
               depth += n;

           }
           return distance(all_reads.reads, max_element(all_reads.reads, all_reads.reads + 4));
        }

        void import_alignment(const BamAlignment& al, const int& pos, const int& bindex){
            uint16_t b_index = base_index(al.AlignedBases[pos]);
            BQ.reads[bindex] += (al.Qualities[pos] - 33);
            MQ.reads[bindex] += (al.MapQuality);
            if(al.IsReverseStrand()){ 
                rev_reads.reads[bindex] += 1; 
            } 
            else{ 
                fwd_reads.reads[bindex] += 1; 
            }
        }
};

 
class ExperimentSiteData{
    public:
        string m_initial_data;
        vector<SampleSiteData> sample_data;
        SampleNames snames;
        char m_ref_base;

        ExperimentSiteData(SampleNames sn, string initial_data, char ref_base){
            snames = sn;
            m_initial_data = initial_data;
            m_ref_base = ref_base;
            //fill constructor doesn't work here?
            for(size_t i= 0; i < sn.size(); i++){
                sample_data.push_back(SampleSiteData());
            }
        }

        void summarize(ostream* out_stream){

            // call genotypes, keep track of each sample and number of each 
            // possible allele
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
                //Looks like messy data. Print out the read matrix so we can
                //unerstand what going on, add an empty line to the output
                *out_stream  << m_initial_data
                   << "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t" << endl;
                cerr << "Skipping " << m_initial_data <<  endl;
                cerr << "Read Matrix:" << endl;
                for (size_t i=0; i < sample_data.size(); i++){
                    cerr << snames[i] << '\t';
                    for (size_t j=0; j<4; j++){
                        cerr << sample_data[i].all_reads.reads[j] << '\t';
                    }
                    cerr << endl;
                }
                return; // 
            }
            auto it = find_if(genotypes.begin(), genotypes.end(), 
                    [&](int v) {return v==mutant_base;});
            uint32_t mutant = distance(genotypes.begin(), it);
            
            //summarise the data from the mutant strain
            SampleSiteData * ms = &sample_data[mutant];
            double f_mutant_m = ms->all_reads.reads[mutant_base]/(double)ms->depth;
            uint16_t* F_m = &ms->fwd_reads.reads[mutant_base];
            uint16_t* R_m = &ms->rev_reads.reads[mutant_base];
            double m_MQs = ms->MQ.reads[mutant_base]/(double)ms->all_reads.reads[mutant_base];
            double m_BQs = ms->BQ.reads[mutant_base]/(double)ms->all_reads.reads[mutant_base];

            int other_MQs_sum = 0;
            int other_BQs_sum = 0;
            int Q_denom = 0;
            for(size_t i = 0; i < 4; ++i){
                if(i != mutant_base){
                    other_MQs_sum += ms->MQ.reads[i];
                    other_BQs_sum += ms->BQ.reads[i];
                    Q_denom += ms->all_reads.reads[i];
                }
            }
            double other_MQs = other_MQs_sum/(double)Q_denom;
            double other_BQs = other_BQs_sum/(double)Q_denom;

            //and now the WT strains
            int N_mutant_wt = 0; //mutant allele freq in strains with WT allel
            int F_wt = 0;// n forward and reverse reads for rference base
            int R_wt = 0;// in wildtype strains
            int wt_MQs_sum = 0;
            int wt_BQs_sum = 0;
            int wt_depth = 0;

            uint16_t ref_bindex  = base_index(m_ref_base);
            for (size_t i=0; i < sample_data.size(); i++){
                if (i != mutant){
                    SampleSiteData* s = &sample_data[i];
                    N_mutant_wt += s->all_reads.reads[mutant_base];

                    F_wt += s->fwd_reads.reads[ref_bindex];
                    R_wt += s->rev_reads.reads[ref_bindex];
    
                    wt_MQs_sum += accumulate(s->MQ.reads, s->MQ.reads + 4, 0);
                    wt_BQs_sum += accumulate(s->BQ.reads, s->BQ.reads + 4, 0);

                    wt_depth += s->depth;

                }
            }
            double wt_MQs = wt_MQs_sum/(double)wt_depth;
            double wt_BQs = wt_BQs_sum/(double)wt_depth;
            double f_mutant_wt = N_mutant_wt/(double)wt_depth;
 
            
            *out_stream  << m_initial_data << '\t'
                         << mutant_base << '\t'
                         << snames[mutant] << '\t'
                         << f_mutant_m << '\t'  //freq. of mutant base in mutant strain
                         << f_mutant_wt << '\t' //freq. mutant base in WTs
                         << *F_m << '\t' << *R_m << '\t' << F_wt << '\t' << R_wt << '\t'
                         << m_BQs << '\t' << other_BQs << '\t' << wt_BQs << '\t'
                         << m_MQs << '\t' << other_MQs << '\t' << wt_MQs << '\t'    
                         << endl;

            }
            
        
};


class FilterVisitor: public PileupVisitor{
    public: 
        FilterVisitor(BamAlignment& ali, 
                      const SamHeader& header,
                      const SampleNames& samples, 
                      ostream* out_stream,
                      int ref_pos,
                      string input_data,
                      char ref_base):
            PileupVisitor(), m_header(header), m_samples(samples), m_ref_pos(ref_pos), 
                             m_out_stream(out_stream), m_initial_data(input_data),
                             m_ref_base(ref_base)
            {  } 
        ~FilterVisitor(void) { }


    public:
        void Visit(const PileupPosition& pileupData){
            if (pileupData.Position != m_ref_pos){
                return;
            }

            auto target_site = ExperimentSiteData(m_samples, m_initial_data, m_ref_base);
            for (auto it =  pileupData.PileupAlignments.begin();
                      it != pileupData.PileupAlignments.end();
                      it++){
                if(it->Alignment.MapQuality > 13){//TODO options for baseQ, mapQ
                    int const *pos = &it->PositionInAlignment;
                    if(it->Alignment.Qualities[*pos] > 46){//TODO user-defined qual cut 
                        uint16_t b_index = base_index(it->Alignment.QueryBases[*pos]);
                        if (b_index < 4){
                            string tag_id;
                            it->Alignment.GetTag("RG", tag_id);
                            string sm = m_header.ReadGroups[tag_id].Sample;
                            uint32_t sindex = find_sample_index(sm,m_samples);
                            target_site.sample_data[sindex].import_alignment(it->Alignment, *pos, b_index);               
                        }
                    }
                }
            }
        target_site.summarize(m_out_stream);
    }

    private:
        SampleNames m_samples;
        ostream* m_out_stream;
        int m_ref_pos;
        char m_ref_base;
        string m_initial_data;
        SamHeader m_header;
       // ExperimentSiteData target_site;
    
};

int main(int argc, char* argv[]){
    string bam_path;
    string input_path;
    vector< string > sample_names;
    namespace po = boost::program_options;
    po::options_description cmd("Command line args");
    cmd.add_options()
        ("help,h", "Print a help message")
        ("bam,b", po::value<string>(&bam_path)->required(), "Path to BAM file")
        ("bam-index,x", po::value<string>()->default_value(""), "Path to BAM index (default is <bam_path>.bai")
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
    
    string index_path = vm["bam-index"].as<string>();
    if(index_path == ""){
        index_path = bam_path + ".bai";
    }

    BamReader experiment;
    experiment.Open(bam_path);
    experiment.OpenIndex(index_path);
    SamHeader header = experiment.GetHeader();
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
        char ref_base = L[i+1];
        int pos = stoul(pos_s);
        int ref_id = experiment.GetReferenceID(chr);
        
        experiment.SetRegion(ref_id, pos, ref_id, pos+1);
        FilterVisitor *f = new FilterVisitor(ali, 
                                             header,
                                             sample_names,
                                             &outfile,
                                             pos, L, ref_base);
        pileup.AddVisitor(f);
        while( experiment.GetNextAlignment(ali) ) {
            pileup.AddAlignment(ali);
        }
        pileup.Flush();
    }
    return 0;
}




            






