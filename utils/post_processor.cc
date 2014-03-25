#include <iostream>
#include <stdint.h>
#include <string>
#include <algorithm>
#include <vector>
#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"

#include "model.h"
#include "parsers.h"

using namespace std;
using namespace BamTools;


class SiteData{

    public:
        string name;
        vector<uint16_t> BQs;
        vector<uint16_t> MQs;
        uint32_t fwd_reads;
        uint32_t rev_reads;
        ReadData base_calls;

        SiteData(string sname){
          name = sname;
          base_calls.key = 0;
        }
                  
        int get_genotype(){
           //These have already been called for mutation-ness, and are haploid
           //so, to make a start, we are just calling the most common base
           return distance(base_calls.reads, max_element(base_calls.reads, base_calls.reads + 4 ));
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


int sum(vector<int> input){
    return(  accumulate(input.begin(), input.end(), 0) );
}

double mean(vector<int> input){
    double s = sum(input);
    return(s/input.size());
}
        

typedef vector<SiteData> SiteDataVector;

class FilterVisitor: public PileupVisitor{
    public: 
        FilterVisitor(BamAlignment ali, ostream *out_stream,
                      SampleNames samples, int ref_pos):

            PileupVisitor(), m_samples(samples), m_ostream(out_stream), 
                             m_ref_pos(ref_pos){
                nsamp = m_samples.size(); 
            } 
        ~FilterVisitor(void) { }
    private:
        SampleNames m_samples;
        int m_ref_pos;
        ostream* m_ostream;
        int nsamp;

    public:
        void Visit(const PileupPosition pileupData){
            SiteDataVector sample_data;
            for (size_t i = 0; i < nsamp; i++){
                sample_data.push_back(SiteData(m_samples[i]));
            }
            for (auto it =  pileupData.PileupAlignments.begin();
                      it != pileupData.PileupAlignments.end();
                      it++){
                if(pileupData.Position == m_ref_pos){
                    int const *pos = &it->PositionInAlignment;
                    if(it->Alignment.Qualities[*pos] > 46){//TODO user-defined qual cut to match first call?
                        uint16_t b_index = base_index(it->Alignment.AlignedBases[*pos]);
                        if (b_index < 4){
                            string tag_id;
                            it->Alignment.GetTag("RG", tag_id);
                            uint32_t sindex = find_sample_index(get_sample(tag_id), m_samples);
                            sample_data[sindex].import_alignment(it->Alignment, *pos, b_index);                        
                        }
                    }
                }
            }

            vector<uint16_t> genotypes;
            int gfreqs[4];
            for(size_t i=0; i < nsamp; i++){
                uint16_t g = sample_data[i].get_genotype();
                gfreqs[g] += 1;
                genotypes.push_back(g);                        
            }
            if(count(gfreqs, gfreqs+4, 1) == 1){//should only be one mutant
                    // OK let's collect all that data....
                    auto it = find_if(genotypes.begin(), genotypes.end(),
                        [](int v) {return v==1;});
                    unsigned mutant = *it;
                    uint16_t mutant_base = genotypes[mutant];
                    int mutant_alleles;
                    int mutant_allele_denom;
                    for (size_t i=0; i < nsamp; i++){
                        if(i != mutant){
                            mutant_alleles += sample_data[i].base_calls.reads[mutant_base];
                            mutant_allele_denom += sample_data[i].fwd_reads + sample_data[i].rev_reads;
                        }
                    }
                    double mutant_freq = (double)mutant_alleles/mutant_allele_denom;
                    *m_ostream << mutant_freq << '\t' 
                               << endl;

            }
        }
    
};
                                 

int main(int argc, char* argv[]){
    return 0;
}
