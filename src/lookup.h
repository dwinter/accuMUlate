/*
 * lookup.h
 *
 *  Created on: 11/25/14
 *      Author: Steven Wu
 */


#pragma once
#ifndef LOOKUP_H_
#define LOOKUP_H_

#include "string"
//#include <pair>
#include <array>


const int ANCESTOR_COUNT = 10; //TODO: Maybe get the constant.h file


namespace LookupTable {
//    class Lookup{
//struct LookupTable {

    extern const int index_vector[4][10]; //Same as summary stat now, number of mismatches between D (o) and A(m,f)
    extern const double summary_stat_index_lookup[4][10];

    extern const double summary_stat_same_lookup_table[10][4];

    extern const int index_converter_16_to_10_single[16];
    extern const int index_converter_4_4_to_10[4][4];

    extern const int index_converter_4_4_to_16[4][4];
    extern const int index_converter_10_to_16[10];
    extern const std::array<int, 2> index_converter_16_to_4_4[16];

    extern const std::string genotype_lookup_10[10];

};
#endif //LOOKUP_H_
