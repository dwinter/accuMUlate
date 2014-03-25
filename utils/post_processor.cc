#include <iostream>
#include <stdint.h>
#include <string>
#include <algorithm>
#include <vector>
#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"

#include "model.h"

using namespace std;

typedef vector<SiteData> SiteDataVector;

class SiteData{

    public:
        string name
        uint16_t genotype
        vector<uint16_t> BQs;
        vector<uint16_t> MQs;
        uint32_t fwd_reads;
        uint32_t rev_reads;
        ReadData base_calls;

        SiteData(string sname){
          name = string;
          base_calls.key = 0;
        }
                  
        void call_genotype(){
           //These have already been called for mutation-ness, and are haploid
           //so, to make a start, we are just calling the most common base
           genotype = distance(base_calls.reads, max_element(base_calls.reads, 4));
        }
        void add_alignment_data(const BamAlignment& al, const int& pos, const int& bindex){
            uint16_t b_index = base_index(al.AlignedBases[*pos]);
            base_calls.reads[bindex] += 1;
            BQs.push_back(al.Qualities[*pos]);
            MQs.push_back(al.MapQuality);
            if(al.IsReverseStrand){ 
                rev_reads += 1; 
            } 
            else{ 
                fwd_reads += 1; 
            }
        }
};



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

class FilterVisitor: public PileupVisitor{
    public: 
        FilterVisitor: public PileupVisitor(BamAlignment ali, 
                                            ostream *out_stream,
                                            SampleNames samples,
                                            int ref_pos)
            PileupVisitor(), m_samples(samples), m_ostream(out_stream), 
                             m_ref_pos(ref_pos){
                nsamp = samples.size(); 
            } 
        ~FilterVisitor(void) { }
    public:
        void Visit(const PileupPosition pileupData){
            SiteDataVector sample_data;
            for (size_t i = 0; i < nsamp i++){
                site.push_back(SiteData(SampleNames[i]))
            }
            for (auto it =  pileupData.PileupAlignments.begin();
                      it != pileupData.PileupAlignments.end();
                      it++){
                if(pileupData.Position == m_ref_pos){
                    int const *pos = &it->PositionInAlignment;
                    if(it->Alignment.Qualities[*pos] > 46){//TODO deal with quality cut
                        uint16_t b_index = base_index(it->Alignment.AlignedBases[*pos]);
                        if (bindex < 4){
                            int sindex = find_sample_index(get_sample(tag_id), m_samples);
                            sample_data[sindex].add_alignment_data(it->Alignment, 
                                                                    *pos, b_index)
                                    
                        }
                    }
                }
            }
            genotypes = vector<uint16_t>;
            for(size_t i=0; i < nsamp; i++){
                genotypes.push_back(sample_data[i].genotype());
            }
            // 
                                        




        

class FreqVisitor: public PileupVisitor{
    public:
        FreqVisitor(BamAlignment ali, string chr, int ref_pos):
            PileupVisitor(), m_ref_pos(ref_pos), m_chr(chr) { }
        ~FreqVisitor(void) { }
    public:
        void Visit(const PileupPosition& pileupData){
            ReadData bcalls;
            bcalls.key = 0;
            for(auto it = begin(pileupData.PileupAlignments);
                     it !=end(pileupData.PileupAlignments);
                     ++it){
                if(pileupData.Position == m_ref_pos){
                    int const *pos = &it->PositionInAlignment;
                    if(it->Alignment.Qualities[*pos] > 46){
                        uint16_t b_index = base_index(it->Alignment.AlignedBases[*pos]);
                        if (b_index < 4){
                            bcalls.reads[b_index] += 1;
                        }
                    }
                }
            }
            
            int m = 0;
            int nreads = 0;
            for(size_t i = 0; i < 4; i++){
                nreads += bcalls.reads[i];
                if(bcalls.reads[i] > m){
                    m = bcalls.reads[i];
                }
            }
            if(nreads > 10){ 
                double result = 1 - (m/(double)nreads);        
                if(result > 0.05){ 
                    cout << m_chr << '\t' << m_ref_pos << '\t' << result << '\t' << nreads << endl;                        
                }
            }
                            
            
        }      
    private:
        int m_ref_pos;
        string m_chr;
};
        

int main(int argc, char* argv[]){
    fstream bed;
    bed.open(argv[1]);
    BamReader bam;
    bam.Open(argv[2]);
    BamAlignment ali;
    PileupEngine pileup;
    string L;
    while(getline(bed, L)){
        string chr_s;
        size_t i = 0;
        for(;L[i] != '\t';i++){
            chr_s.push_back(L[i]);
        }
        int chr = bam.GetReferenceID(chr_s);


        i+= 1;
        string pos_s;
        for(;L[i] != '\t'; i++){
            pos_s.push_back(L[i]);
        }
        int pos = stoul(pos_s);
        bam.SetRegion(chr,pos, chr, pos);
        FreqVisitor *f = new FreqVisitor(ali, chr_s, pos);
        pileup.AddVisitor(f); 
        while(bam.GetNextAlignment(ali)){
            pileup.AddAlignment(ali);
        }
    };
    pileup.Flush();
    return 0;
}
