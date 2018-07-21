#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <limits>

#include <api/BamReader.h>
#include "model.h"
#include "parsers.h"


using namespace std;
using namespace BamTools;


//Base visitor, which does two things per site:
//Set a flag descrbing
ReadDataVisitor::ReadDataVisitor(LocalBamToolsUtils::Fasta &idx_ref,
                                 SampleMap &samples,
                                 int qual_cut, int mapping_cut) :
        PileupVisitor(),
        m_idx_ref(idx_ref), m_samples(samples), m_mapping_cut(mapping_cut) {

    qual_cut_char = (char) (qual_cut + 33);

    uint32_t max = 0;
    for (auto item : m_samples) {
        if (item.second != MAX_UINT32 && item.second > max) {
            max = item.second;
        }
    }
    total_sample_count = max + 1; //Plus ref sindex==0;

    ReadDataVector bcalls (total_sample_count, ReadData{0});
    site_data =  {0, bcalls};
}


bool ReadDataVisitor::GatherReadData(const LocalBamToolsUtils::PileupPosition &pileupData) {

    int pos  = pileupData.Position;
    // Pileup produces all alignments that _overlap_ with a region. If we have
    // set a region for variant calling then we will also get some of the
    // reads upstream and downstream of that region. So we need to check a given
    // position is actually in our region-of-interest
    if( pos < region_start || pos >= region_end){  
        return false;
    }
    //If we are in the right region, start processing the data
    m_idx_ref.GetBase(pileupData.RefId, pos, current_base);
    uint16_t ref_base_idx = base_index_lookup[(int) current_base];
    if( ref_base_idx > 4 ) { // TODO: This treats all non IUPAC codes as masks. Document this is we keep it
        return false;
    }

    site_data.reference = ref_base_idx;
    ReadDataVector &bcalls = site_data.all_reads;
    std::fill(bcalls.begin(), bcalls.end(), ReadData(0));

    for (auto it = begin(pileupData.PileupAlignments); it != end(pileupData.PileupAlignments); ++it) {
        
        if(not it->IsCurrentDeletion ){
            int32_t pos_in_alignment = it->PositionInAlignment;
            if (include_site(it->Alignment, pos_in_alignment, m_mapping_cut, qual_cut_char)) {
                uint32_t sindex = GetSampleIndex(it->Alignment);
    
                if (sindex != MAX_UINT32) {
                    uint16_t bindex = base_index_lookup[(int) it->Alignment.QueryBases[pos_in_alignment]];
                    if (bindex < 4) {
                        bcalls[sindex].reads[bindex] += 1;
                    }
                }
            }
        }
    }

    return true;
};


SiteStatsSummary ReadDataVisitor::CalculateStats(const LocalBamToolsUtils::PileupPosition &pileupData, int mutant_index, int mallele_index, int rotate_by){

    SiteStatsSummary res = {};
    if(mallele_index == -1){
        // It is possible to get "mutations" that make no sense (A->A) or are
        // not easy to handle (GT->CA), skip these cases and document this
        // behaviour.
        return res;
    }
    int pos  = pileupData.Position;
    m_idx_ref.GetBase(pileupData.RefId, pos, current_base);
    uint16_t ref_base_idx = base_index_lookup[(int) current_base];

    //Collect the site's data
    SiteStatsData raw_data = {};
    for (auto it = begin(pileupData.PileupAlignments); it != end(pileupData.PileupAlignments); ++it) {
        if(not it->IsCurrentDeletion ){
            int32_t pos_in_alignment = it->PositionInAlignment;
            if (include_site(it->Alignment, pos_in_alignment, m_mapping_cut, qual_cut_char)) {            
                uint32_t sindex = GetSampleIndex(it->Alignment);
                uint16_t bindex = base_index_lookup[(int) it->Alignment.QueryBases[pos_in_alignment]];
                if ( sindex == mutant_index + 1) {
                    bindex += rotate_by;
                    if(bindex > 3){
                        bindex -= 4;
                    }
                }
                if (sindex != MAX_UINT32) {
                    bool rev = it->Alignment.IsReverseStrand();
                    bool paired = raw_data.Paired_anc += it->Alignment.IsMateMapped();
                    if (mallele_index == bindex){
                        raw_data.Paired_mutant += paired;
                        raw_data.MQ_mutant.push_back(it->Alignment.MapQuality);
                        if (paired){
                            raw_data.Insert_mutant.push_back(it->Alignment.InsertSize);
                        }                
                        if (sindex == mutant_index + 1){//Mutant index is among _descendants_
                            if (rev){                   // GetSampleIndex includes the ancestor as sample zero
                                raw_data.MM_R += 1;                            
                            }
                            else {
                                raw_data.MM_F += 1; 
                            }
                        }
                        else{
                            raw_data.AM += 1;
                        }
                    }
                    else{
                        raw_data.Paired_anc += paired;
                        raw_data.MQ_anc.push_back(it->Alignment.MapQuality);
                        if (paired){
                            raw_data.Insert_anc.push_back(it->Alignment.InsertSize);
                        }                
                        if (sindex == mutant_index + 1){                                           
                            if (rev){
                                raw_data.MO_R += 1;
                            }
                            else{
                                raw_data.MO_F += 1;
                            }
                        }
                        else{
                            raw_data.AA += 1;
                        }
                    }
                        
                }        
            }
        }
    }
    //Do the stats
    int N_minor_mutant = raw_data.MO_R + raw_data.MO_F;
    int unpaired_anc = (N_minor_mutant + raw_data.AA) - raw_data.Paired_anc;
    int unpaired_mutant = (raw_data.MM_F + raw_data.MM_R + raw_data.AM) - raw_data.Paired_mutant;
    res.MQ_AD = ad_two_sample_test(raw_data.MQ_mutant, raw_data.MQ_anc);
    res.Insert_AD = ad_two_sample_test(raw_data.Insert_mutant, raw_data.Insert_anc);
    res.FisherPairBias = fisher_exact_test(raw_data.Paired_mutant, unpaired_mutant, raw_data.Paired_anc, unpaired_anc);
    res.FisherStrandBias = fisher_exact_test(raw_data.MM_R, raw_data.MM_F, raw_data.MO_R, raw_data.MO_F );
    res.DP= raw_data.MM_R + raw_data.MM_F + raw_data.AM + raw_data.AA + N_minor_mutant;
    res.NM_F = raw_data.MM_F;
    res.NM_R = raw_data.MM_R;
    res.NM_WT = raw_data.AM;
    res.N_minor = N_minor_mutant;
    return(res);

}

uint32_t ReadDataVisitor::GetSampleIndex(const BamAlignment &alignment) {
    std::string tag_id_2;
    if (! alignment.GetTag("RG", tag_id_2)) {
        std::cout << "ERROR: Could not find RG tag in tag data\t" << alignment.TagData << std::endl;
    }
    return m_samples[tag_id_2];
}

void ReadDataVisitor::SetRegion(BedInterval target_region){
    region_start = target_region.start;
    region_end = target_region.end;    
}


BedFile::BedFile(string bed_file_name) {
    bed_file.open(bed_file_name);
}

int BedFile::get_interval(BedInterval &current_interval) {
    string L;
    if (getline(bed_file, L)) {
        stringstream Lstream(L);
        string chrom;
        string start_s;
        string end_s;
        getline(Lstream, chrom, '\t');
        getline(Lstream, start_s, '\t');
        getline(Lstream, end_s, '\t');
        current_interval = BedInterval{chrom,
                                       stoul(start_s),
                                       stoul(end_s)};
        return 0;
    }
    else { return 1; }

}





bool include_site(const BamAlignment &alignment, const int &pos, const uint16_t &map_cut, const char &qual_cut) {
//    const BamAlignment *ali = &(pileup.Alignment);
    if (alignment.MapQuality > map_cut) {
        char reference = alignment.Qualities[pos];    
        if (reference > qual_cut) {           
                    //Is primary line (not supp. or secondary ali)
            return ( (alignment.AlignmentFlag & 0x900) == 0  && not (alignment.IsDuplicate()) && not (alignment.IsFailedQC()) );
        }
    }

    return false;
}



const int base_index_lookup[128] ={
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
//	        A  B  C  D  e  f  G  H  i  j  K  l  M  N  o
        -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,    // 64-79
//   	 p  q  R  S  T  U  V  W  x  Y  z
        -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,	// 80-95
//	        A  B  C  D  e  f  G  H  i  j  K  l  M  N  o
        -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,   	// 96-111
//  	 p  q  R  S  T  U  V  W  x  Y  z
        -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1		// 112-127
};


