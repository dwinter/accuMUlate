/*
 * DirichletMultinomialDistribution.h
 *
 *  Created on: Nov 7, 2014
 *      Author: Steven Wu
 */

#ifndef DIRICHLETMULTINOMIALDISTRIBUTION_H_
#define DIRICHLETMULTINOMIALDISTRIBUTION_H_


#include <data_struct.h>

class DirichletMultinomialDistribution {
public:
	DirichletMultinomialDistribution();
	virtual ~DirichletMultinomialDistribution();

//	double DirichletMultinomialLogProbability(double alphas[4], ReadData data);
//	void DirichletMultinomialRandom(double alhpas[4], ReadData &data);


};

double DirichletMultinomialLogProbability(double (&alphas)[4], ReadData const &data);

#endif /* DIRICHLETMULTINOMIALDISTRIBUTION_H_ */
