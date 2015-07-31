#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <limits>       

#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"
#include "utils/bamtools_fasta.h"


#include "model.h"
#include "parsers.h"


using namespace std;
using namespace BamTools;



//Base visitor, which does two things per site:
//Set a flag descrbing
ReadDataVisitor::ReadDataVisitor(
                        const RefVector& bam_references, 
                        Fasta& idx_ref,
                        SampleMap& samples, 
                        const ModelParams& p,  
                        BamAlignment& ali, 
                        int qual_cut,
                        int mapping_cut):

            PileupVisitor(),
                m_bam_references(bam_references),            
                m_idx_ref(idx_ref),
                m_samples(samples), 
                m_params(p),
                m_ali(ali), 
                m_qual_cut(qual_cut),  
                m_mapping_cut(mapping_cut) {}
//        ~ReadDataVisitor(void) { }
//    public:

bool ReadDataVisitor::GatherReadData(const PileupPosition& pileupData) {
    //Like it says, collect a sites reads. If the site is good to call
    //from set the site_data object and return `true`.  
    uint64_t pos  = pileupData.Position;
    m_idx_ref.GetBase(pileupData.RefId, pos, current_base);
    uint16_t ref_base_idx = base_index(current_base);
    if( ref_base_idx > 4 ) { // TODO: This treats all non IUPAC codes as masks. Document this is we keep it 
        return false;
    }
    ReadDataVector bcalls (m_samples.size(), ReadData{{ 0,0,0,0 }}); 
    for(auto it = begin(pileupData.PileupAlignments);
             it !=  end(pileupData.PileupAlignments); 
             ++it){
        if( include_site(*it, m_mapping_cut, m_qual_cut) ){
            it->Alignment.GetTag("RG", tag_id);
            uint32_t sindex = m_samples[tag_id]; 
                if( sindex  != std::numeric_limits<uint32_t>::max()  ){
                uint16_t bindex  = base_index(it->Alignment.QueryBases[it->PositionInAlignment]);
                if (bindex < 4 ){
                    bcalls[sindex].reads[bindex] += 1;
                }
            }
        }
    }
    site_data =  {ref_base_idx, bcalls};
    return true;
};


BedFile::BedFile(string bed_file_name){
     bed_file.open(bed_file_name);
}

int BedFile::get_interval(BedInterval& current_interval){
    string L;
    if(getline(bed_file, L)){
        stringstream Lstream(L);
        string chrom;
        string start_s;
        string end_s;        
        getline(Lstream, chrom, '\t');
        getline(Lstream, start_s, '\t');
        getline(Lstream, end_s, '\t');
        current_interval = BedInterval{ chrom, 
                                         stoul(start_s), 
                                         stoul(end_s)};
        return 0;
    }
   else {return 1;}
    
}





//
//Helper functions for parsing data out of BAMs

uint16_t base_index(char b){
    switch(b){        
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
             return -1 ;
         default: // Unknown base, alert in debug mode?
             return -1;
    }
}


//uint32_t find_sample_index(string s, SampleNames sv){
//    for (size_t i=0; i < sv.size(); i++){
//        if(s.compare(sv[i])==0){
//            return i;
//        }
//    }
//    return(-1); //TODO refactor this to  update sample in place
//}

bool include_site(PileupAlignment pileup, uint16_t map_cut, uint16_t qual_cut){
    const BamAlignment *ali = &pileup.Alignment;
    if(ali->MapQuality > map_cut){
        uint16_t bqual = static_cast<short>(ali->Qualities[pileup.PositionInAlignment]) - 33;
        if(bqual > qual_cut){
            return(not (ali->IsDuplicate()) && not(ali->IsFailedQC()) && ali->IsPrimaryAlignment());
        }
    }
    return false;
}


//int main(){
//    return 0;
//}


//int main() {
//    FastaReference reference_g ("test/test.fai");
//    int chr_idx;
//    reference_g.get_ref_id("scf_8254727", chr_idx);
//    cout << "chr index = " << chr_idx << " (should be 9)" << endl;
//
//    BedFile bed ("test/test.bed");
//    BedInterval current_line;
//    while(bed.get_interval(current_line) == 0){
//         cout << current_line.chr << endl;
//    }
//    return 0;
//}
//
   
