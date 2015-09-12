/*
 * DirichletMultinomialDistribution.cpp
 *
 *  Created on: Nov 7, 2014
 *      Author: Steven Wu
 */

#include "DirichletMultinomialDistribution.h"

//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_rng.h>
#include <cmath>
#include <iostream>


//#include <cmath>
//#include <iostream>
//#include "../model.h"

DirichletMultinomialDistribution::DirichletMultinomialDistribution() {
	// TODO Auto-generated constructor stub

}

DirichletMultinomialDistribution::~DirichletMultinomialDistribution() {
	// TODO Auto-generated destructor stub
}


//
//double DirichletMultinomialDistribution::DirichletMultinomialLogProbability(double alphas[4], ReadData data) {
//	// TODO: Cache most of the math here
//	// TODO: Does not include the multinomail coefficient
//
////	for (int i = 0; i < 4; ++i) {
////			data.reads[i] = 0;
////			cout << alphas[i] << endl;
////			alphas[i] = 0.25;
////		}
////	data.reads[0] = 3;
//
//	int read_count = data.reads[0]+data.reads[1]+data.reads[2]+data.reads[3];
//
////	printf("RC:%d\n",read_count);
//	double alpha_total = alphas[0]+alphas[1]+alphas[2]+alphas[3];
//	double result = 0.0;
//	for(int i : {0,1,2,3}) {
//		for(int x = 0; x < data.reads[i]; ++x) {
//			result += log(alphas[i]+x);
////			cout << alphas[i]+x << "\t";
////			cout << lgamma(alphas[i]+x) << endl;
//
//		}
//		std::cout<< result << "\n";
//	}
////	double sum = 0;
//	for(int x = 0; x < read_count; ++x){
//		result -= log(alpha_total+x);
////		sum += log(alpha_total+x);
//	}
//	std::cout<< result << "\n";
////	exit(-1);
////	printf("sum denom:%f\n",sum);
//	return result;
//}

//void DirichletMultinomialDistribution::DirichletMultinomialRandom(
//		double alphas[4], ReadData& data) {
//
//	  const gsl_rng_type * T;
//	  gsl_rng * r;
//	  gsl_rng_env_setup();
//
//	  T = gsl_rng_mt19937;
//	  r = gsl_rng_alloc (T);
//
//	  int i, n = 10;
//	  double mu = 3.0;
//
//	  /* create a generator chosen by the
//	     environment variable GSL_RNG_TYPE */
//
//	  double theta[4] = {0.1,0.1,0.1,0.1};
////	gsl_rng_
//	gsl_ran_dirichlet(r, 10, alphas, theta);
//	/*
//	simPop
//function (J = 10, K = 20, n, pi, theta)
//{
//    if (length(n) == 1)
//        n <- rep(n, J)
//    if (missing(pi))
//        pi <- rnorm(K, mean = 14, sd = 4)
//    else K <- length(pi)
//    pi <- pi/sum(pi)
//    P <- rdirichlet(J, pi * (1 - theta)/theta)
//    X <- matrix(0, J, K)
//    for (i in 1:J) X[i, ] <- rmultinom(1, n[i], P[i, ])
//    list(theta = theta, pi = pi, data = X)
//}

//>
//
//	*/
//
//
//}




double DirichletMultinomialLogProbability(double (&alphas)[4], ReadData const &data) {
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
	for(int x = 0; x < read_count; ++x){
		result -= log(alpha_total+x);
	}

	return result;
}
