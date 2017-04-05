
#include <iostream>
#include <cstdint>
#include <cmath>
#include <vector>
#include <initializer_list>
#include <iostream>
#include <map>
#include <fstream>
#include <memory>
#include <iomanip>
#include "Eigen/Dense"
#include "Eigen/StdVector"
#include "model.h"

using namespace std;



//Create the mutation matrices which represent the probabilites of different
//histories with or without mutation. 
TransitionMatrix F81(const ModelParams &params) {
	double beta = 1.0;
	for(auto d : params.nuc_freq)
		beta -= d*d;
	beta = 1.0/beta;
	beta = exp(-beta*params.mutation_rate);

	TransitionMatrix m;
	for (int i : { 0, 1, 2, 3 }) {
		double prob = params.nuc_freq[i]*(1.0-beta);
		m.row(i) = Eigen::Vector4d::Constant(prob);
	}
	m.diagonal() += Eigen::Vector4d::Constant (beta);

	return m;
}

MutationMatrix MutationAccumulation(const ModelParams &params, bool and_mut){
	TransitionMatrix m = F81(params);
    if(params.ploidy_ancestor == 1){//haploid->haploid design
        if(!and_mut){
            return m;
        }
        MutationMatrix result = MutationMatrix(4,4);           
        //identity matrix initialzer is for Matrices only, these are (despite teh
        //name) arrays so buid it up:
        for(int i : {0,1,2,3}){
            for (int j : {0,1,2,3}){
                if(i != j){
                    result(i,j) = m(i,j);
                }
                else{
                    result(i,j) = 0.0;
                }
            }
        }
        return result;
    }
    if(params.ploidy_descendant == 1){//diploid -> haploid design
        MutationMatrix result = MutationMatrix(16,4);
        for(int i : {0,1,2,3}) {
		    for(int j : {0,1,2,3}) {
			    for(int k : {0,1,2,3}) {
				    result(i*4+j,k) = 0.0;
        				if(!and_mut || i != k)
					result(i*4+j,k) += 0.5*m(i,k);
		        		if(!and_mut || j != k)
					result(i*4+j,k) += 0.5*m(j,k);
                }
            }
	    }
	    return result;
    }
    //only diploid->diploid left 
    MutationMatrix result = MutationMatrix(16,16);
	for(int i : {0,1,2,3}) {
		for(int j : {0,1,2,3}) {
			for(int k : {0,1,2,3}) {
    			for(int l : {0,1,2,3}) {                
    				result(i*4+j,k*4+l) = 0.0;
    				if(!and_mut || i != k || j != l ){//TODO: check this is right transition prob
					    result(i*4+j,k*4+l) += m(i,k) * m(j,l);
                    }
                }
            }
        }
    }
    return result;
}


GenotypeVector diploid_genotypes = { 
    {"AA", 1, 0}, 
    {"AC", 2, 1}, 
    {"AG", 2, 2}, 
    {"AT", 2, 3}, 
    {"CC", 1, 5}, 
    {"CG", 2, 6}, 
    {"CT", 2, 7}, 
    {"GG", 1, 10}, 
    {"GT", 2, 11}, 
    {"TT", 1, 15}, 
};



GenotypeVector haploid_genotypes = { 
    {"A", 1, 0}, 
    {"C", 1, 1}, 
    {"G", 1, 2}, 
    {"T", 1, 3}, 
};


//Find mutant allele from genotype->genotype transmissions. Ignore apparent
//double mutations (-1).
static int diploid_mutant_matrix[16][16] = {
    // Starred columns are already covered by an earlier row, present here
    // to keep 'mat-index" from GenotypeProperties object
    // 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
    //AA  AC  AG  AT  CA* CC  CG  CT  GA* GC* GG  GT  TA* TC* TG* TT 
    {-1,  1,  2,  3, -1,  1, -1, -1, -1, -1,  2, -1, -1 ,1, -1,  3}, //AA.
    { 0, -1,  2,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1 ,1, -1, -1}, //AC.      
    { 0,  1, -1,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1 ,1, -1, -1}, //AG.
    { 0,  1,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 ,1, -1, -1}, //AT.
    { 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 ,1, -1, -1}, //CA*
    { 0,  0, -1, -1, -1, -1,  2,  3, -1, -1,  2, -1, -1 ,1, -1,  3}, //CC.
    {-1,  0,  0, -1, -1,  1, -1,  3, -1, -1, -1, -1, -1 ,1, -1, -1}, //CG.
    {-1,  0, -1,  0, -1,  1,  2, -1, -1, -1, -1, -1, -1 ,1, -1,  3}, //CT.
    {-1, -1, -1,  3, -1, -1,  2, -1, -1, -1,  2,  3, -1 ,1, -1, -1}, //GA*
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 ,1, -1, -1}, //GC*
    { 0, -1,  0, -1, -1,  1,  1, -1, -1, -1, -1,  3, -1 ,1, -1,  3}, //GG
    {-1, -1,  0,  0, -1, -1,  3,  2, -1, -1,  2, -1, -1 ,1, -1,  3}, //GT
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 ,1, -1, -1}, //TA*
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 ,1, -1, -1}, //TC*
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 ,1, -1, -1}, //TG*
    { 0, -1, -1,  0, -1,  1, -1,  1, -1, -1,  2,  2, -1 ,1, -1, -1}, //TT
};





GenotypeProbs PopulationProbs(SequencingFactory &sf, int ref_allele, int ploidy_ancestor) {
	if(ploidy_ancestor==2){
		DiploidProbs result = sf.getRefDiploidProbs(ref_allele);
		return result;
	}
	HaploidProbs result = sf.getRefHaploidProbs();
	return result;
}

GenotypeProbs Sequencing(SequencingFactory &sf, ReadData data, int ploidy) {
	if(ploidy == 2){
		DiploidProbs result = sf.GetDiploidSequencing(data);
		return result;
	}
	HaploidProbs result = sf.GetHaploidSequencing(data);
	return result;
}

double TetMAProbability(const ModelParams &params, SequencingFactory &sf,
						const ModelInput &site_data,
						const MutationMatrix &m, const MutationMatrix &mn) {


	auto it = site_data.all_reads.begin();
	GenotypeProbs pop_genotypes = PopulationProbs(sf, site_data.reference, params.ploidy_ancestor);
	GenotypeProbs anc_genotypes = Sequencing(sf, *it, params.ploidy_ancestor);

	anc_genotypes *= pop_genotypes;
	GenotypeProbs num_genotypes = anc_genotypes;

	for(++it; it != site_data.all_reads.end(); ++it) {
		GenotypeProbs p = Sequencing(sf, *it, params.ploidy_descendant);

		anc_genotypes *= (m.matrix()*p.matrix()).array();
		num_genotypes *= (mn.matrix()*p.matrix()).array();
	}

	return 1.0 - num_genotypes.sum()/anc_genotypes.sum();
}

double TetMAProbOneMutation(const ModelParams &params, SequencingFactory &sf,
							const ModelInput &site_data,
							const MutationMatrix &m, const MutationMatrix &mn) {

	auto it = site_data.all_reads.begin();

	GenotypeProbs pop_genotypes = PopulationProbs(sf, site_data.reference, params.ploidy_ancestor);
	GenotypeProbs anc_genotypes = Sequencing(sf, *it, params.ploidy_ancestor);
	anc_genotypes *= pop_genotypes;
	GenotypeProbs denom = anc_genotypes;   //product of p(Ri|A)
	GenotypeProbs nomut_genotypes = anc_genotypes; //Product of p(Ri & noMutatoin|A)
	GenotypeProbs mut_genotypes = anc_genotypes;      //Sum of p(Ri&Mutation|A=x)
	mut_genotypes.setZero();

	for(++it; it != site_data.all_reads.end(); ++it) {
		GenotypeProbs p = Sequencing(sf, *it, params.ploidy_descendant);
		GenotypeProbs dgen =  (mn.matrix()*p.matrix()).array();
		GenotypeProbs agen = (m.matrix()*p.matrix()).array();
		nomut_genotypes *= dgen;
		mut_genotypes += (agen/dgen - 1); //(agen+dgen)/agen
		denom *= agen;
	}
	double result = (nomut_genotypes * mut_genotypes).sum() / denom.sum();
	return(result);
}


MutationDescription DescribeMutant(const ModelParams &params, SequencingFactory &sf, const ModelInput site_data, const MutationMatrix m, const MutationMatrix mn) {
    MutationMatrix mt = m - mn;
    int ndesc = site_data.all_reads.size() - 1 ;
    GenotypeVector from_genotypes = (params.ploidy_ancestor == 2 ? diploid_genotypes : haploid_genotypes);
    GenotypeVector to_genotypes = (params.ploidy_descendant == 2 ? diploid_genotypes : haploid_genotypes);

    std::vector<GenotypeProbs, Eigen::aligned_allocator<Eigen::ArrayXd> > lower_mn(ndesc);
    std::vector<GenotypeProbs, Eigen::aligned_allocator<Eigen::ArrayXd> > lower_m(ndesc);
	GenotypeProbs pop_genotypes = PopulationProbs(sf, site_data.reference, params.ploidy_ancestor);	
    GenotypeProbs anc_genotypes = Sequencing(sf, site_data.all_reads[0], params.ploidy_ancestor);
	anc_genotypes *= pop_genotypes;
    GenotypeProbs denom = anc_genotypes;
    double max_mu = 0;
    uint16_t mutant_line;
    //Calculate P(G|R), store as a matrix
	for(size_t i = 1 ; i <= ndesc; ++i) {
        GenotypeProbs p = Sequencing(sf, site_data.all_reads[i], params.ploidy_descendant);
        GenotypeProbs agen = (m.matrix() * p.matrix()).array();
        lower_mn[i-1] = (mn.matrix() * p.matrix()).array();
        lower_m[i-1] = agen;
        denom *= agen;
    }
    //Find the line with the highest probabilty of being the _only_ mutant
	for(size_t i = 0 ; i < ndesc; ++i) {
        GenotypeProbs mut = anc_genotypes;
	    for(size_t j = 0 ; j < ndesc; ++j) {
            if( i == j){
                mut *= lower_m[j]  - lower_mn[j];
            } else {
                mut *= lower_mn[j];
            }
        }
        double p_one_mutation=  mut.sum() /denom.sum();
        if(p_one_mutation > max_mu){
            max_mu = p_one_mutation;
            mutant_line = i; // Get back to sample-name index
        }
    }
    //For that line, what is the most likely genotype change.
    GenotypeProbs mutant_sequencing =  Sequencing(sf, site_data.all_reads[mutant_line+1], params.ploidy_descendant);
    string from, to;
    double mu = 0;
    double line_denom = 0;
    uint16_t mutant_allele = 0;
    for( GenotypeProperties A : from_genotypes) {
        for( GenotypeProperties D : to_genotypes) {
            double res = denom[ A.mat_index ] * A.ways * mutant_sequencing[ D.mat_index ] * D.ways;
            line_denom += res;
            if(res > mu){
                from = A.bases;
                to = D.bases;
                mu = res;
                if(params.ploidy_descendant == 1){
                    mutant_allele = D.mat_index;
                }
                else{
                    mutant_allele = diploid_mutant_matrix[A.mat_index][D.mat_index];

                }
            }
        }
    }
    MutationDescription final = {mutant_line, mutant_allele, from, to, max_mu, mu/line_denom , denom.sum()};
    return final;
     
}

