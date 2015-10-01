//
// Created by steven on 4/3/15.
//

#include "unittest_utils.h"

const double ERROR_THRESHOLD = 1e-14;//std::numeric_limits<double>::epsilon()*10;


void ASSERT_GENOTYPES(HaploidProbs expected, HaploidProbs data) {
    for (int i = 0; i < expected.size(); ++i) {
        ASSERT_NEAR(expected[i], data[i], ERROR_THRESHOLD );
    }
}

void ASSERT_GENOTYPES(DiploidProbs expected, DiploidProbs data) {
    for (int i = 0; i < expected.size(); ++i) {
        ASSERT_NEAR(expected[i], data[i], ERROR_THRESHOLD );
    }
}

void ASSERT_GENOTYPES(DiploidProbsIndex10 expected, DiploidProbsIndex10 data) {
    for (int i = 0; i < expected.size(); ++i) {
        ASSERT_NEAR(expected[i], data[i], ERROR_THRESHOLD );
    }
}
//TODO: should be a  template for this??
//template <typename T>
//void ASSERT_VECTOR_LIKE(T expected, T data) {
//    ASSERT_EQ(expected.size(), data.size() );
//    for (int i = 0; i < expected.size(); ++i) {
//        ASSERT_NEAR(expected[i], data[i], ERROR_THRESHOLD );
//    }
//}
