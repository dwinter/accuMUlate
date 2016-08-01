#ifndef stats_H
#define stats_H

#include <algorithm>
#include <iostream>
#include <cmath>
#include <vector>


double fisher_exact_test(int a11, int a12, int a21, int a22);   

double ad_two_sample_test(std::vector<int> a, std::vector<int> b);

double phred(double p);



#endif
