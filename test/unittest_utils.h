//
// Created by steven on 4/3/15.
//

#ifndef _UNITTEST_UTILS_H_
#define _UNITTEST_UTILS_H_

#include "model.h"
#include "gtest/gtest.h"
//#include "constant.h"



extern const double ERROR_THRESHOLD;

extern void ASSERT_GENOTYPES(HaploidProbs expected, HaploidProbs data);

extern void ASSERT_GENOTYPES(DiploidProbs expected, DiploidProbs data);

extern void ASSERT_GENOTYPES(DiploidProbsIndex10 expected, DiploidProbsIndex10 data);


template <typename T>
void ASSERT_VECTOR_LIKE(T expected, T data) {
    ASSERT_EQ(expected.size(), data.size() );
    for (int i = 0; i < expected.size(); ++i) {
        ASSERT_NEAR(expected[i], data[i], ERROR_THRESHOLD );
    }
}


#endif //__UNITTEST_UTILS_H_
