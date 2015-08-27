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


bool ReadDataVisitor::GatherReadDataNew(const PileupPosition& pileupData) {

    //Like it says, collect a sites reads. If the site is good to call
    //from set the site_data object and return `true`.
    uint64_t pos  = pileupData.Position;

    m_idx_ref.GetBase(pileupData.RefId, pos, current_base);
    uint16_t ref_base_idx = base_index2[(int) current_base];
    if( ref_base_idx > 4 ) { // TODO: This treats all non IUPAC codes as masks. Document this is we keep it
        return false;
    }
    ReadDataVector bcalls (m_samples.size(), ReadData{{ 0,0,0,0 }});
    for(auto it = begin(pileupData.PileupAlignments);
        it !=  end(pileupData.PileupAlignments); ++it){
//        if( include_site(*it, m_mapping_cut, m_qual_cut) ){

        int32_t pos_in_alignment = it->PositionInAlignment;
        if (include_site_4(it->Alignment, pos_in_alignment, m_mapping_cut, qual_cut_char)) {
//            it->Alignment.GetTag("RG", tag_id);
//            uint32_t sindex = m_samples[tag_id];

            int sindex = GetSampleIndex(it->Alignment.TagData);

            if( sindex  != std::numeric_limits<uint32_t>::max()  ){
//                uint16_t bindex  = base_index(it->Alignment.QueryBases[it->PositionInAlignment]);
                uint16_t bindex = base_index2[(int) it->Alignment.QueryBases[pos_in_alignment]];
                if (bindex < 4 ){
                    bcalls[sindex].reads[bindex] += 1;
                }
            }
        }
    }
    site_data =  {ref_base_idx, bcalls};
    return true;
};


int ReadDataVisitor::GetSampleIndex(std::string const &tag_data) {
    size_t start_index = tag_data.find(rg_tag);
    if (start_index != std::string::npos) {
        start_index += 4;
    }
    else {
        size_t x = tag_data.find("RG");//TODO: Check for (char)0, RG, Z
        std::cout << "ERRER: " << tag_data << "\t" << x << std::endl;
    }
    size_t end_index = tag_data.find(ZERO_CHAR, start_index);
    std::string tag_id_2 = tag_data.substr(start_index, (end_index - start_index));
//
//
////            string tag_id;
////            it->Alignment.GetTag("RG", tag_id);
////            string sm = m_header.ReadGroups[tag_id].Sample;
////            if(tag_id != tag_id_2){
////                cout << tag_id << "\t" << tag_id_2 <<  "\t" << start_index << "\t" << end_index << "\n" << tag_data << endl;
////                exit(-1);
////            }
////            if(sm != sm2){
////                cout << sm << "\t" << sm2 << endl;
////                exit(-2);
////            }
//
////            string sm2 = map_tag_sample_two_stage[tag_id_2]; //TODO:catch exception
////            uint32_t sindex = m_samples[sm2]; //TODO check samples existed!
//
    int sindex = m_samples[tag_id_2];
    return sindex;
}


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


bool include_site_4(const BamAlignment & alignment, const int &pos, const uint16_t &map_cut, const char &qual_cut){
//    const BamAlignment *ali = &(pileup.Alignment);
    if(alignment.MapQuality > map_cut){
        char reference = alignment.Qualities[pos];

        if(reference > qual_cut){
            return(not (alignment.IsDuplicate()) && not(alignment.IsFailedQC()) && alignment.IsPrimaryAlignment());
        }
    }

    return false;
}



int base_index2[128] ={
//    17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,	// 0-15
//    17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,	// 16-31
////                                          -
//    17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,	// 32-47
////                                                ?
//    17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,16,	// 48-63
////	    A  B  C  D  e  f  G  H  i  j  K  l  M  N  o
//    17, 0,11, 1,12,16,16, 2,13,16,16,10,16, 7,15,16,	// 64-79
////	 p  q  R  S  T  U  V  W  x  Y  z
//    16,16, 5, 9, 3, 3,14, 8,16, 6,16,17,17,17,17,17,	// 80-95
////	    A  B  C  D  e  f  G  H  i  j  K  l  M  N  o
//    17, 0,11, 1,12,16,16, 2,13,16,16,10,16, 7,15,16,	// 96-111
////	 p  q  R  S  T  U  V  W  x  Y  z
//    16,16, 5, 9, 3, 3,14, 8,16, 6,16,17,17,17,17,17		// 112-127

        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,	// 0-15
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,	// 16-31
//                                          -
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,	// 32-47
//                                                ?
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,	// 48-63
//	    A  B  C  D  e  f  G  H  i  j  K  l  M  N  o
        -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,    // 64-79
//	 p  q  R  S  T  U  V  W  x  Y  z
        -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,	// 80-95
//	    A  B  C  D  e  f  G  H  i  j  K  l  M  N  o
        -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,   	// 96-111
//	 p  q  R  S  T  U  V  W  x  Y  z
        -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1		// 112-127
};



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
   

