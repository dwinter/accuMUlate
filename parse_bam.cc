#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"
#include <iostream>
#include <stdint.h>
#include <map>
#include <vector>
#include <string>

#include "model.h"

using namespace std;
using namespace BamTools;

string bases = "ACGT";

typedef vector< string > SampleNames;


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
    return(13); //TODO refactor this to  update sampe in place!
}

size_t base_index(char b){
    switch(b){
        case 'A':
            return 0;
         case 'T':
             return 3;
         case 'C':
             return 1;
         case 'G':
             return 2;
         case '-':
         case 'N':
             return string::npos;
         default:
             cerr << "Don't know what to make of " << b <<endl;
             return string::npos;
    };
}

                             
class VariantVisitor : public PileupVisitor{
    public:
        VariantVisitor(const RefVector& references, const SampleNames samples, ModelParams p, int q=13): 
            PileupVisitor(), m_references(references), m_samples(samples), m_q(q), m_params(p)
            {
                nsamp = m_samples.size();
            }
        ~VariantVisitor(void) { }
    public:
         void Visit(const PileupPosition& pileupData) {
             cout << m_references[pileupData.RefId].RefName << "\t";
             cout << pileupData.Position << "\t";
             cout << m_references << '\t';
             ReadDataVector bcalls (nsamp, ReadData{{ 0,0,0,0 }}); //fill constructor
             string tag_id;
             for(auto it = begin(pileupData.PileupAlignments);
                      it != end(pileupData.PileupAlignments); 
                      ++it){
                 int const *pos = &it->PositionInAlignment;
                 if (it->Alignment.Qualities[*pos] - 33 > m_q){
                     it->Alignment.GetTag("RG", tag_id);
                     int sindex = find_sample_index(get_sample(tag_id), m_samples);
                     size_t bindex  = base_index(it->Alignment.AlignedBases[*pos]);
                     if (bindex != string::npos){
                         bcalls[sindex].reads[bindex] += 1;
                     }
                 }
            }
            ModelInput d = {"", 1, 0, bcalls};
            cout << TetMAProbOneMutation(m_params,d) << '\t'  
                 << TetMAProbability(m_params,d) << endl;          
        }
    private:
        RefVector m_references;
        SampleNames m_samples;
        int m_q;
        int nsamp;
        ModelParams m_params;
        
            
};



int main(){
    BamReader myBam; 
    myBam.Open("scf_8254670.bam");
    RefVector references = myBam.GetReferenceData();
    SampleNames all_samples {"M0", "M19", "M20", "M28","M25", "M29", 
                             "M40", "M44","M47", "M50","M51", "M531"};
     ModelParams p = { 
       0.001, 
       {0.25, 0.25, 0.25, 0.25}, 
       1.0e-8,
       0.01,
       0.001,
       0.001,
    };
    BamAlignment ali;
    PileupEngine pileup;
    VariantVisitor *v = new VariantVisitor(references, all_samples,p);
    pileup.AddVisitor(v);
    while( myBam.GetNextAlignment(ali)){
        pileup.AddAlignment(ali);
        ` 
    };
    pileup.Flush();
    return 0;
}
        
        



