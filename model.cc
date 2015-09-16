
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


//int main(){
//      ModelParams p = { 
//        0.0001, 
//        {0.38, 0.12, 0.12, 0.38}, 
//        1e-2,
//        0.01,
//        0.01,
//        0.05, 
//        1,1
//};
//   MutationMatrix mt = MutationAccumulation(p, true);
//   MutationMatrix m = MutationAccumulation(p, false);
//   cerr << m << endl << endl;
//   cout << m -mt << endl;
//   return 0;
////}
//// Uncommon and compile with this:
//// clang++ -std=c++11 -Ithird-party/bamtools/src/ -Lboost_progam_options model.cc
//
// to play around with / debug results.
//int main(){
//    ModelParams p = { 
//        0.0001, 
//        {0.38, 0.12, 0.12, 0.38}, 
//        1e-8,
//        0.01,
//        0.01,
//        0.05,
//        2,2
//    };
//   MutationMatrix mt = MutationAccumulation(p, true);
//   MutationMatrix m = MutationAccumulation(p, false);
//   MutationMatrix mn = m - mt;
//    ModelInput two_vars = { 2, 
//        {
//        { 0, 30,  0,  0},
//        { 0, 30,  0,  0},
//        { 0, 30,  0,  0},
//        { 0, 30,  0,  0},
//        { 0,  0,  0, 30},
//        { 0, 0,   0, 30},
//        { 0, 30,  0,  0},
//        { 0, 30,  0,  0},
//        { 0, 30,  0,  0},
//        { 0, 30,  0,  0},
//        { 0, 30,  0,  0},
//        { 0, 30,  0,  0},
//
//        }
//    };
//
//
//    ModelInput  one_vars = { 1,
//        {
//        { 0, 30,  0,  0},
//        { 0,  0,  0, 30},
//        { 0, 30,  0,  0},
//        { 0, 30,  0,  0},
//        { 0, 30,  0,  0},
//        { 0, 30,  0, 0}}
//    };    
//
//   ModelInput no_vars = {1,
//      {{0,   0,  0, 5},       
//       {3,   0,  0, 4},       
//       {8,   0,  0, 0},
//       {3,   0,  0, 3},     
//       {0,  0,  0, 2},       
//       {2,  0,  0, 3},
//       {0,  0,  0, 0},  
//       {0,  0,  0, 5},       
//       {0,  0,  0, 4},
//       {1,  0,  0, 3},
//       {0,  0,  0, 28},
//       {0 ,  0,  0, 3}}
//   };
//    
//    cout << "___With the no-variant data___" << endl;
////    cout << "P(one|data)= "<<  TetMAProbOneMutation(p,no_vars)<< endl;
////    cout << "P(any|data)= " << TetMAProbability(p,no_vars) << endl;
//
//    cout << TetMAProbability(p,no_vars,m,mn) << endl;
//    
//    cout << "___With the one-variant data___" << endl;  
////    cout << "P(one|data)= "<<  TetMAProbOneMutation(p,one_vars)<< endl;
////    cout << "P(any|data)= " << TetMAProbability(p,one_vars) << endl;
//      cout << TetMAProbability(p,one_vars,m,mn) << endl;
//    
//    cout << "___With the two-variant data___" << endl;
////    cout << "P(one|data)= "<<  TetMAProbOneMutation(p,two_vars)<< endl;
//    cout << "P(any|data)= " << TetMAProbability(p,two_vars,m,mn) << endl;   
//
//    
////    cout << "calculating the same number once: " << TetMAProbOneMutation(p,two_vars) << endl;
////    cout << "then another time: " << TetMAProbOneMutation(p,two_vars) << endl;
////    TetMAProbOneMutation(p,no_vars);
////    cout << "And once more after calling from the the no-vars data: " << TetMAProbOneMutation(p,two_vars) << endl;
//    return 0;
//}
//






