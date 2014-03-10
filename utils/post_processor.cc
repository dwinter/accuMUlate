#include <iostream>
#include <stdint.h>
#include <map>
#include <vector>
#include <string>
//#include <basic_string>
#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"

#include "model.h"

using namespace BamTools;
using namespace std;

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
        FreqVisitor(BamAlignment ali):
            PileupVisitor()  { }
        ~FreqVisitor(void) { }
    public:
        void Visit(const PileupPosition& pileupData){
            ReadData bcalls;
            bcalls.key = 0;
            for(auto it = begin(pileupData.PileupAlignments);
                     it !=end(pileupData.PileupAlignments);
                     ++it){
                int const *pos = &it->PositionInAlignment;
                if(it->Alignment.Qualities[*pos] > 46){
                    uint16_t b_index = base_index(it->Alignment.AlignedBases[*pos]);
                    if (b_index < 4){
                        bcalls.reads[b_index] += 1;
                    }
                }
            }
            for(size_t i = 0; i < 4; i++){
                cout << bcalls.reads[i] << '\t';
            }
            cout << endl;
        }
    private:
};
        

int main(int argc, char* argv[]){
    fstream bed;
    bed.open(argv[1]);
    BamReader bam;
    bam.Open(argv[2]);
    BamAlignment ali;
    PileupEngine pileup;
    FreqVisitor *f = new FreqVisitor(ali);
    pileup.AddVisitor(f);
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
        bam.SetRegion(chr,pos, chr, pos+1);
        while(bam.GetNextAlignment(ali)){
            pileup.AddAlignment(ali);
        };
    };
    pileup.Flush();
    return 0;
}
