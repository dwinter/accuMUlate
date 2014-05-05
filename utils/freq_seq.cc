#include <iostream>
#include <stdint.h>
#include <string>
#include <algorithm>

#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"

#include "model.h"
#include "parsers.h"

using namespace BamTools;
using namespace std;


class FreqVisitor: public PileupVisitor{
    public:
        FreqVisitor(BamAlignment ali, 
                    const SamHeader& header,
                    const SampleNames& samples,
                    int ref_pos,
                    char ref_base,
                    string initial_data
                    ):
            PileupVisitor(), m_ref_pos(ref_pos), m_samples(samples), 
                             m_header(header), m_ref_base(ref_base),
                             m_initial_data(initial_data) {
                nsamp = m_samples.size(); 
            }
        ~FreqVisitor(void) { }
    public:
        void Visit(const PileupPosition& pileupData){
            uint16_t ref_base_index = base_index(m_ref_base); 
            if(pileupData.Position != m_ref_pos){
                return;
            } 

            vector<int> non_ref(nsamp, 0);
            vector<int> denoms(nsamp, 0);
            for(auto it = begin(pileupData.PileupAlignments);
                     it !=end(pileupData.PileupAlignments);
                     ++it){
                if(it->Alignment.MapQuality > 13){
                    int const *pos = &it->PositionInAlignment;
                    if(it->Alignment.Qualities[*pos] > 46){///TODO user-define cut-offs
                        uint16_t b_index = base_index(it->Alignment.QueryBases[*pos]);
                        if (b_index < 4){
                            string tag_id;
                            it->Alignment.GetTag("RG", tag_id);
                            string sm = m_header.ReadGroups[tag_id].Sample;
                            uint32_t sindex = find_sample_index(sm, m_samples);
                            if(b_index != ref_base_index){
                                non_ref[sindex] += 1;
                            }
                            denoms[sindex] += 1;
                        }
                    }   
                }
            }
            cout << m_initial_data << '\t';
            for(size_t i = 0; i < nsamp; i++){// will return -nan for sample with no coverage.
                cout << non_ref[i]/(double)denoms[i] << '\t';
            }
            cout << endl;
            return; // don't nead the rest of the
        } 
    private:
        int nsamp;
        SampleNames m_samples;
        char m_ref_base;
        SamHeader m_header;
        int m_ref_pos;
        string m_initial_data;
};
        
        

int main(int argc, char* argv[]){
    fstream bed;
    bed.open(argv[1]);
    string bam_path = argv[2];
    BamReader bam;
    bam.Open(bam_path);
    bam.OpenIndex(bam_path + ".bai");
    SamHeader header = bam.GetHeader();
    SampleNames samples;
    for(auto it = header.ReadGroups.Begin(); it!= header.ReadGroups.End(); it++){
        if(it->HasSample()){
                samples.push_back(it->Sample);
        }
    }
    BamAlignment ali;
    string L;
    while(getline(bed, L)){
        PileupEngine pileup;
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
        char ref_base = L[i+1];

        bam.SetRegion(chr,pos, chr, pos+1);

        FreqVisitor *f = new FreqVisitor(ali, 
                                         header,
                                         samples,
                                         pos,
                                         ref_base,
                                         L);
        pileup.AddVisitor(f); 
        while(bam.GetNextAlignment(ali)){
            pileup.AddAlignment(ali);
        }
    pileup.Flush();
    }
    return 0;
}
      
