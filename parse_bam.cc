#include <iostream>
#include <stdint.h>
#include <map>
#include <vector>
#include <string>

#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"
#include "utils/bamtools_fasta.h"

#include "model.h"

using namespace std;
using namespace BamTools;

typedef vector< string > SampleNames;
typedef map < string, int> ChrMap;

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

//int get_id_by_name(string chr, const Fasta& refseq){
//    int result = find_if(begin(refseq.FastaPrivate.Index), 
//                         end(refseq.FastaPrivate.Index),
//                         [&](const FastaIndexData  &fid) {fid.Name == chr});
//   return result;
//}


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

                             
class VariantVisitor : public PileupVisitor{
    public:
        VariantVisitor(const RefVector& bam_references, 
                       const Fasta& fasta_references,
                       SampleNames samples, 
                       const ModelParams& p,  
                       BamAlignment& ali, 
                       int q =13,
                       double prob_cutoff=0.1):

            PileupVisitor(), m_idx_ref(fasta_references), m_bam_ref(bam_references), 
                             m_samples(samples), m_q(q), m_params(p), m_ali(ali)
            {
                nsamp = m_samples.size();
            }
        ~VariantVisitor(void) { }
    public:
         void Visit(const PileupPosition& pileupData) {
             string chr = m_bam_ref[pileupData.RefId].RefName;
             uint64_t pos  = pileupData.Position;
//             int chr_idx = get_id_by_name(pileupData.RefId, m_idx_ref);
             m_idx_ref.GetBase(pileupData.RefId, pos, current_base);
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
                     if (bindex < 4){
                         bcalls[sindex].reads[bindex] += 1;
                     }
                 }
            }
            ModelInput d = {"", 1, base_index(current_base), bcalls};
            double prob = TetMAProbability(m_params,d);
            if(prob >= 0.0){//NOTE: Probablity cut off is hard-coded atm
             cout << chr << '\t' << pos << '\t' << current_base << '\t' << 
                 prob << '\t' << TetMAProbOneMutation(m_params,d) << endl;          
            }
        }
    private:
        Fasta m_idx_ref; 
        RefVector m_bam_ref;
        SampleNames m_samples;
        int m_q;
        int nsamp;
        ModelParams m_params;
        BamAlignment& m_ali;
        char current_base;
        uint64_t chr_index;
        
            
};



int main(){
    BamReader myBam; 
    myBam.Open("test/test.bam");
    RefVector references = myBam.GetReferenceData(); 
    cerr << "buliding fasta index..." << endl;
    Fasta reference_genome;
    reference_genome.Open("test/test.fasta");
    reference_genome.CreateIndex("test/test.fai");
 //   ifstream idx ("test/test.fai");

    //bamtools api makes the actual index private
    //only public acess to a sequence is via it ID (index)
    //so need to make a map of them
//    string idx_line;
//    ChrMap idx_map;
//
//    while(getline(idx, idx_line)){
//        size_t i = 0;
//        string name;
//        int counter = 0;
//        for(; L[i] != '\t', i++){
//            name.push_back(L[i])
//        }
//        idx_map[name] = counter;
//        counter += 1;
//    }
//


            
    reference_genome.Close();
    reference_genome.Open("test/test.fasta", "test/test.fai");

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
    VariantVisitor *v = new VariantVisitor(references, reference_genome, all_samples,p, ali);
    pileup.AddVisitor(v);
    cerr << "calling variants" << endl;
    while( myBam.GetNextAlignment(ali)){
        pileup.AddAlignment(ali);
         
    };
    pileup.Flush();
    return 0;
}
        
        



