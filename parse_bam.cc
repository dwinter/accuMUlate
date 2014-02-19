#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"
#include <iostream>
#include <stdint.h>
#include <map>
#include <vector>
#include <string>

using namespace std;
using namespace BamTools;

union ReadData {
    uint16_t reads[4];
    uint64_t key;
};





string bases = "ACGT";

string get_sample(string& tag){
    string res;
    for(size_t i = 0; tag[i] != '_'; i++) {
        res.push_back(tag[i]);
    }
    return(res);            
}




class VariantVisitor : public PileupVisitor{
    public:
        VariantVisitor(const RefVector& references): 
            PileupVisitor(), m_references(references)
            { }
        ~VariantVisitor(void) { }
    public:
         void Visit(const PileupPosition& pileupData){
             cout << m_references[pileupData.RefId].RefName << "\t";
             cout << pileupData.Position << "\t";
             cout << pileupData.PileupAlignments.size() << "\t";
             map<string, ReadData> sample_counts {
                {"M0",   ReadData{{0,0,0,0}}},
                {"M19",  ReadData{{0,0,0,0}}},
                {"M20",  ReadData{{0,0,0,0}}},
                {"M28",  ReadData{{0,0,0,0}}},
                {"M40",  ReadData{{0,0,0,0}}},
                {"M44",  ReadData{{0,0,0,0}}}, 
                {"M50",  ReadData{{0,0,0,0}}},
                {"M531", ReadData{{0,0,0,0}}},
             };
             string tag_id;   
             for(auto it = begin(pileupData.PileupAlignments);
                      it != end(pileupData.PileupAlignments); 
                      ++it){
                   
                 it->Alignment.GetTag("RG", tag_id);
                 string s = get_sample(tag_id);
                 char b = it->Alignment.AlignedBases[it->PositionInAlignment];
                 sample_counts[s].reads[bases.find(b)] += 1;        
                 if(pileupData.Position == 98){
                     cout << s  ;
                     cout << b << " ";
                 }
             }
             cout << endl;                
        }
    private:
        RefVector m_references;
};





int main(){
    BamReader myBam; 
    myBam.Open("scf_8254670.bam");
    RefVector references = myBam.GetReferenceData();
    cout << references.size() << endl;
    cout << references[1].RefLength<< endl;
    BamAlignment ali;
    PileupEngine pileup;
    VariantVisitor *v = new VariantVisitor(references);
    pileup.AddVisitor(v);
    while( myBam.GetNextAlignment(ali)){
        pileup.AddAlignment(ali); 
    };
    pileup.Flush();
    return 0;
}
        
        



