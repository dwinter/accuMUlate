#include <iostream>
#include <stdint.h>
#include <map>
#include <vector>
#include <string>
#include <algorithm>

#include "boost/program_options.hpp"

#include "boost_input_utils.h"
#include "model.h"
#include "denom.h"

using namespace std;
using namespace BamTools;



//bool include_sample(const ModelParams &params, const ReadDataVector fwd, const ReadDataVector rev,  const ReadDataVector site_data, int sindex, uint16_t ref_base, double pcut, bool central){
//    
//    //Can't be included if you don't have 3fwd, 3rev so check that before we do
//    //any number crunching    
//    //if(central == true){
//        auto Fm = max_element(fwd[sindex].reads, fwd[sindex].reads+4);
//        if (*Fm < 3){
//            return false;
//        }
//        
//        auto Rm = max_element(rev[sindex].reads, rev[sindex].reads+4);
//        if (*Rm < 3){
//            return false;
//        }
//   // }
//    ReadDataVector _site_data = site_data;
//    for(size_t i = 1; i < 4; i++){
//        rotate( begin(_site_data[sindex].reads), begin(_site_data[sindex].reads) + i, end(_site_data[sindex].reads) );
//        ModelInput d = { ref_base, _site_data };
//        double p  = TetMAProbability(params, d); 
//        if (p > 0.1){
//            return true;
//        }
//    }
//    return false;
//}
//    


DenomVisitor::DenomVisitor( const RefVector& bam_references, 
                            const SamHeader& header,
                            LocalBamToolsUtils::Fasta& idx_ref,
                            SampleMap& samples,
                            BamAlignment& ali, 
                            int qual_cut,
                            int mapping_cut,
                            DenomCounterVector &denoms,
                            ModelParams& params,
                            DenomCriteria criteria,                            
                            double prob_cut):

    ReadDataVisitor(idx_ref, samples, qual_cut, mapping_cut),
            
        m_bam_ref(bam_references), 
        m_header(header),m_prob_cut(prob_cut),sf(m_params), 
        m_ali(ali), m_denoms(denoms), m_params(params), m_criteria(criteria) {
           
           
        MutationMatrix mt = MutationAccumulation(m_params, true);
        m_mut_paths = MutationAccumulation(m_params, false);
        m_nomut_paths = m_mut_paths - mt;
        nsamp = samples.size(); 
    }



bool DenomVisitor::MeetsCriteria(SiteStatsSummary stats){
    if(stats.DP < m_criteria.DP_min){
        return(false);
    }
    if(stats.DP > m_criteria.DP_max){
        return(false);
    }
    if(stats.NM_F < m_criteria.mutant_strand_min){
        return(false);
    }
    if(stats.NM_R < m_criteria.mutant_strand_min){
        return(false);
    }
    if(stats.N_minor > m_criteria.ancestral_in_mutant_max){
        return(false);
    }
    if(stats.NM_WT > m_criteria.mutant_in_ancestral_max){
        return(false);
    }
    if(stats.MQ_AD > m_criteria.MQ_AD_max){
        return(false);
    }
    if(stats.Insert_AD > m_criteria.Insert_AD_max){
        return(false);
    }
    if(stats.FisherStrandBias < m_criteria.Fisher_std_min){
        return(false);
    }
    if(stats.FisherPairBias < m_criteria.Fisher_map_min){
        return(false);
    }
    return true;
}


         //For every site update the m_denoms by incrementing the entry for 
         //each sample meeting the criteria 
void DenomVisitor::Visit(const LocalBamToolsUtils::PileupPosition &pileupData) {
    if (GatherReadData(pileupData)) {
        for(size_t i  = 1; i < nsamp; i++){
            for(size_t b = 0; b < 4; b++){
                ModelInput _site_data = site_data;
                std::rotate( std::begin(_site_data.all_reads[i].reads), std::begin(_site_data.all_reads[i].reads) + b, std::end(_site_data.all_reads[i].reads) );
                double p  = TetMAProbability(m_params,sf, _site_data,  m_mut_paths, m_nomut_paths); 
                if( p > m_prob_cut){
                    MutationDescription details = DescribeMutant(m_params, sf, site_data, m_mut_paths, m_nomut_paths);
                    SiteStatsSummary stats = CalculateStats(pileupData, details.mutant_line, details.mutant_allele_index);
                    if( MeetsCriteria(stats) ){
                        m_denoms[i].counts[_site_data.reference] += 1;
                        break;
                    }
                }
            }
        }
    }
}
 



int main(int argc, char** argv){


    boost::program_options::variables_map vm;
    BoostUtils::ParseDenominateCommandline(argc, argv,vm);
    
    BamReader experiment;                                                                                               
    RefVector references;                                                                                
    SamHeader header;
    LocalBamToolsUtils::Fasta reference_genome;
    
    ModelParams params = BoostUtils::CreateModelParams(vm);    
    BoostUtils::ExtractInputVariables(vm, experiment, references, header, reference_genome);                                                  
    SampleMap samples = BoostUtils::ParseSamples(vm, header);


    DenomCriteria criteria = {
        vm["min-depth"].as<uint32_t>(),
        vm["max-depth"].as<uint32_t>(),
        vm["min-mutant-strand"].as<uint32_t>(),
        vm["max-anc-in-mutant"].as<uint32_t>(),
        vm["max-mutant-in-anc"].as<uint32_t>(),
        vm["max-MQ-AD"].as<double>(),
        vm["max-insert-AD"].as<double>(),
        vm["min-strand-pval"].as<double>(),
        vm["min-mapping-pval"].as<double>()
    };
                                                                                                                                                  
    LocalBamToolsUtils::PileupEngine pileup;                                                                                                  
    BamAlignment ali;
    DenomCounterVector denoms (samples.size(), {0,0,0,0} );



    DenomVisitor *v = new DenomVisitor(
            references,
            header,
            reference_genome,
            samples,
            ali,
            vm["qual"].as<int>(),
            vm["mapping-qual"].as<int>(),
            denoms,
            params,
            criteria,
            vm["prob"].as<double>()
        );
    pileup.AddVisitor(v);
   
    if (vm.count("intervals")){                                                                                                               
        BedFile bed (vm["intervals"].as<string>());                                                                                           
        BedInterval region;                                                                                                                   
        while(bed.get_interval(region) == 0){                                                                                                 
            int ref_id = experiment.GetReferenceID(region.chr);                                                                               
            experiment.SetRegion(ref_id, region.start, ref_id, region.end);                                                                   
            v->SetRegion(region);                                                                                                             
            while( experiment.GetNextAlignment(ali) ){                                                                                        
                pileup.AddAlignment(ali);                                                                                                     
            }                                                                                                                                 
        pileup.Flush();                                                                                                                       
        }                                                                                                                                     
    }                                                                                                                                         
    else{                                                                                                                                     
        while( experiment.GetNextAlignment(ali)){                                                                                             
            pileup.AddAlignment(ali);                                                                                                         
        }                                                                                                                                     
    }                                                                                                                                         
    pileup.Flush(); 
    
    for(size_t i = 1; i < samples.size(); i++){
        for( size_t j = 0; j < 4; j++){
            cout << denoms[i].counts[j] << '\t'; 
        }
    }
    cout << endl;
    return 0;
}


