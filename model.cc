#include <iostream>
#include <cstdint>
#include <cmath>
#include <vector>
#include <initializer_list>
#include <iostream>
#include <map>
#include <fstream>
#include <memory>
#include "Eigen/Dense"

#include "model.h"


using namespace std;

double DirichletMultinomialLogProbability(double alphas[4], ReadData data) {
	// TODO: Cache most of the math here
	// TODO: Does not include the multinomail coefficient
	int read_count = data.reads[0]+data.reads[1]+data.reads[2]+data.reads[3];
	double alpha_total = alphas[0]+alphas[1]+alphas[2]+alphas[3];
	double result = 0.0;
	for(int i : {0,1,2,3}) {
		for(int x = 0; x < data.reads[i]; ++x) {
			result += log(alphas[i]+x);
		}
	}
	for(int x = 0; x < read_count; ++x)
		result -= log(alpha_total+x);
	return result;
}

//'Population' refers to the pop. from which the ancestor comes
// We deal only with GC/AT bias and (in the diploid case) heterozygosity


DiploidProbs DiploidPopulation(const ModelParams &params, int ref_allele) {
	ReadData d;
	DiploidProbs result;
	double alphas[4];
	for(int i : {0,1,2,3})
		alphas[i] = params.theta*params.nuc_freq[i];
	for(int i : {0,1,2,3}) {
		for(int j=0;j<i;++j) {
			d.key = 0;
			d.reads[ref_allele] = 1;
			d.reads[i] += 1;
			d.reads[j] += 1;
			result[i+j*4] = DirichletMultinomialLogProbability(alphas, d);
			result[j+i*4] = result[i+j*4];
		}
		d.key = 0;
		d.reads[ref_allele] = 1;
		d.reads[i] += 2;
		result[i+i*4] = DirichletMultinomialLogProbability(alphas, d);
	}
	return result.exp();
}

GenotypeProbs PopulationProbs(const ModelParams &params, int ref_allele){
    if(params.ploidy_ancestor==2){
        DiploidProbs result = DiploidPopulation(params, ref_allele);
        return result;
    }
    HaploidProbs result;
    for( int i :{0,1,2,3}) {
       result[i] = params.nuc_freq[i];
    }
    return result;
}


//Create the mutation matrices which represent the probabilites of different
//histories with or without mutation. 

TransitionMatrix F81(const ModelParams &params){
	double beta = 1.0;
	for(auto d : params.nuc_freq)
		beta -= d*d;
	beta = 1.0/beta;
	beta = exp(-beta*params.mutation_rate);
	TransitionMatrix m;
	for(int i : {0,1,2,3}) {
		for(int j : {0,1,2,3}) {
			m(i,j) = params.nuc_freq[i]*(1.0-beta);
		}
		m(i,i) += beta;
    }
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

//Calculate the genotype likelihoods. Again we have one function with returns 
//differently sized results depending on the ploidy of experiment.


DiploidProbs DiploidSequencing(const ModelParams &params, int ref_allele, ReadData data) {
	DiploidProbs result;
	double alphas_total = (1.0-params.phi_diploid)/params.phi_diploid;
	for(int i : {0,1,2,3}) {
		for(int j=0;j<i;++j) {
			double alphas[4];
			for(int k : {0,1,2,3}) {
				if(k == i || k == j)
					alphas[k] = (0.5-params.error_prob/3.0)*alphas_total;
				else
					alphas[k] = (params.error_prob/3.0)*alphas_total;
			}
			result[i*4+j] = DirichletMultinomialLogProbability(alphas, data);
			result[j*4+i] = result[i*4+j];
		}
		double alphas[4];
		for(int k : {0,1,2,3}) {
			if(k == i)
				alphas[k] = (1.0-params.error_prob)*alphas_total;
			else
				alphas[k] = params.error_prob/3.0*alphas_total;
		}
		result[i*4+i] = DirichletMultinomialLogProbability(alphas, data);
	}
	double scale = result.maxCoeff();
	return (result - scale).exp();
}

HaploidProbs HaploidSequencing(const ModelParams &params, int ref_allele, ReadData data) {
	HaploidProbs result;
	double alphas_total = (1.0-params.phi_haploid)/params.phi_haploid;
	for(int i : {0,1,2,3}) {
		double alphas[4];
		for(int k : {0,1,2,3}) {
			if(k == i)
				alphas[k] = (1.0-params.error_prob)*alphas_total;
			else
				alphas[k] = params.error_prob/3.0*alphas_total;
		}
		result[i] = DirichletMultinomialLogProbability(alphas, data);
	}
	double scale = result.maxCoeff();
	return (result - scale).exp();
}

GenotypeProbs Sequencing(const ModelParams &params, int ref_allele, ReadData data, int ploidy) {
    if(ploidy == 2){
        DiploidProbs result = DiploidSequencing(params, ref_allele, data);
        return result;
    }
    HaploidProbs result = HaploidSequencing(params, ref_allele, data);
    return result;
}


double TetMAProbability(const ModelParams &params, const ModelInput site_data, const MutationMatrix m, const MutationMatrix mn) {
    GenotypeProbs pop_genotypes = PopulationProbs(params, site_data.reference);
	auto it = site_data.all_reads.begin();
    GenotypeProbs anc_genotypes = Sequencing(params, site_data.reference, *it, params.ploidy_ancestor);
	anc_genotypes *= pop_genotypes;
    GenotypeProbs num_genotypes = anc_genotypes;
	for(++it; it != site_data.all_reads.end(); ++it) {
		GenotypeProbs p = Sequencing(params, site_data.reference, *it, params.ploidy_descendant);
		anc_genotypes *= (m.matrix()*p.matrix()).array();
		num_genotypes *= (mn.matrix()*p.matrix()).array();
	}

	return 1.0 - num_genotypes.sum()/anc_genotypes.sum();
}

double TetMAProbOneMutation(const ModelParams &params, const ModelInput site_data, const MutationMatrix m, const MutationMatrix mn) {
	GenotypeProbs pop_genotypes = PopulationProbs(params, site_data.reference);	
	auto it = site_data.all_reads.begin();
    GenotypeProbs anc_genotypes = Sequencing(params, site_data.reference, *it, params.ploidy_ancestor);
	anc_genotypes *= pop_genotypes;
  	GenotypeProbs denom = anc_genotypes;   //product of p(Ri|A)
    GenotypeProbs nomut_genotypes = anc_genotypes; //Product of p(Ri & noMutatoin|A)
    GenotypeProbs mut_genotypes = anc_genotypes;      //Sum of p(Ri&Mutation|A=x)
    mut_genotypes.setZero();
	for(++it; it != site_data.all_reads.end(); ++it) {
        GenotypeProbs p = Sequencing(params, site_data.reference, *it, params.ploidy_descendant);
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






