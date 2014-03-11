#include <iostream>
#include <stdint.h>
#include <string>
#include <algorithm>
#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"

#include "model.h"

using namespace BamTools;
using namespace std;

//TODO this is copy-pasted from parse_bam.cc. If is used more 
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
