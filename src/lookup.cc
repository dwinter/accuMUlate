/*
 * lookup.cc
 *
 *  Created on: 11/25/14
 *      Author: Steven Wu
 */

#include "lookup.h"

namespace LookupTable {
//struct LookupTable {
//    class Lookup {
        const int index_vector[4][10] = { //Same as summary stat now, number of mismatches between D (o) and A(m,f)
                {0, 1, 1, 1, 2, 2, 2, 2, 2, 2},// A
                {2, 1, 2, 2, 0, 1, 1, 2, 2, 2},// C
                {2, 2, 1, 2, 2, 1, 2, 0, 1, 2},// G
                {2, 2, 2, 1, 2, 2, 1, 2, 1, 0} // T
        };
        const double summary_stat_index_lookup[4][10] = {
                {0, 1, 1, 1, 2, 2, 2, 2, 2, 2},// A
                {2, 1, 2, 2, 0, 1, 1, 2, 2, 2},// C
                {2, 2, 1, 2, 2, 1, 2, 0, 1, 2},// G
                {2, 2, 2, 1, 2, 2, 1, 2, 1, 0} // T

//        {0, 0.5, 0.5, 0.5, 1, 1, 1, 1, 1, 1},// A
//        {1, 0.5, 1, 1, 0, 0.5, 0.5, 1, 1, 1},// C
//        {1, 1, 0.5, 1, 1, 0.5, 1, 0, 0.5, 1},// G
//        {1, 1, 1, 0.5, 1, 1, 0.5, 1, 0.5, 0} // T
        };

        const double summary_stat_same_lookup_table[10][4] = {
                //A    C    G    T
                {1.0, 0, 0, 0},    //	AA
                {0.5, 0.5, 0, 0},    //	AC
                {0.5, 0, 0.5, 0},    //	AG
                {0.5, 0, 0, 0.5},  //	AT
                {0, 1.0, 0, 0},    //	CC
                {0, 0.5, 0.5, 0},    //	CG
                {0, 0.5, 0, 0.5},  //	CT
                {0, 0, 1.0, 0},    //	GG
                {0, 0, 0.5, 0.5},  //	GT
                {0, 0, 0, 1.0},    //	TT
        };

        const int index_converter_16_to_10_single[16] = {
                0, 1, 2, 3,
                1, 4, 5, 6,
                2, 5, 7, 8,
                3, 6, 8, 9,
        };
        const int index_converter_4_4_to_10[4][4] = {
                {0, 1, 2, 3},
                {1, 4, 5, 6},
                {2, 5, 7, 8},
                {3, 6, 8, 9}
        };
        const int index_converter_4_4_to_16[4][4] = {
                {0, 4, 8,  12},
                {1, 5, 9,  13},
                {2, 6, 10, 14},
                {3, 7, 11, 15}
        };
        const int index_converter_10_to_16[10] = {
                0, 1, 2, 3,
                5, 6, 7,
                10, 11,
                15
        };
        const std::string genotype_lookup_10[10] = {
                "AA", "AC", "AG", "AT",
                "CC", "CG", "CT",
                "GG", "GT", "TT"
        };

        const std::array<int, 2> index_converter_16_to_4_4[16] = {
                {{0,0}},{{1,0}},{{2,0}},{{3,0}},
                {{0,1}},{{1,1}},{{2,1}},{{3,1}},
                {{0,2}},{{1,2}},{{2,2}},{{3,2}},
                {{0,3}},{{1,3}},{{2,3}},{{3,3}}
        };




};

