#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <sstream>

#include "parsers.h"


using namespace std;
using namespace BamTools;

FastaReference::FastaReference(string ref_file_name){
    ifstream ref_file (ref_file_name);
    string L;
    uint64_t cummulative_len = 0;
    while( getline(ref_file, L)){
        size_t i = 0;
        string chrom;
        for(; L[i] != '\t'; i++){
            chrom.push_back(L[i]);
        }
        i += 1;
        string s_len;
        for (; L[i] !='\t'; i++){
            s_len.push_back(L[i]);
        }
        uint32_t chrom_len = stoul(s_len);
        cummulative_len += chrom_len;
        chromosomes.push_back(FastaReferenceData{ chrom, 
                                                  chrom_len, 
                                                  cummulative_len
                                                 });
    }
}
    
void FastaReference::get_ref_id(string search_name, int& chr_id){
    auto it = find_if(chromosomes.begin(), chromosomes.end(),
                     [&search_name](const FastaReferenceData& chrom){
                        return chrom.name == search_name;
                      });
    int idx  = distance(chromosomes.begin(), it);
    chr_id = idx;
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





   
