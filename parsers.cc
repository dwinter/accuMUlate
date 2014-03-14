#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <sstream>

#include "parsers.h"

using namespace std;

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
    auto it = find_if(chromosomes.begin(), 
                     chromosomes.end(),
                     [&search_name](const FastaReferenceData& chrom){
                        return chrom.name == search_name;
                      });
    int idx  = distance(chromosomes.begin(), it);
    chr_id = idx;
}


BedFile::BedFile(string bed_file_name){
    ifstream bed_file (bed_file_name);
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
        getline(Lstream, start_s, '\t');
        current_interval = BedInterval{ chrom, 
                                         stoul(start_s), 
                                         stoul(end_s)};
        return 0;
    }
    else {return 1;}
}



int main() {return 0;}






   
