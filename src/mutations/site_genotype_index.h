/*
 * site_genotypes_index.cc
 *
 *  Created on: 4/5/15
 *      Author: Steven Wu
 */

//
// Created by steven on 4/5/15.
//
#pragma once
#ifndef _SITE_GENOTYPES_INDEX_H_
#define _SITE_GENOTYPES_INDEX_H_

#include <vector>
#include <iostream>
#include <algorithm>
//#include <mutations/model.h>
//#include "mutation_prob.h"
//#include "site_prob.h"
//#include "site_genotypes.h"
//#include "distributions/DirichletMultinomialDistribution.h"
//#include "evolution_models/EvolutionModel.h"
//#include "lookup.h"
//#include "sequence_prob_v1.h"
//#include "constant.h"

class SiteGenotypesIndex {
public:
//    SiteGenotypesIndex (SiteGenotypesIndex&& other) = default;
    SiteGenotypesIndex (SiteGenotypesIndex&& other) :ancestor_index(other.ancestor_index) {
//        std::cout << "move constructor" << std::endl;
        other.ancestor_index = 0;
        descendant_index = std::move(other.descendant_index);

    }
    SiteGenotypesIndex& operator= (SiteGenotypesIndex&& other) {
//        std::cout << "move assignment" << std::endl;
        ancestor_index = other.ancestor_index;
        other.ancestor_index = 0;
        descendant_index = std::move(other.descendant_index);

        return *this;
    }

//    SiteGenotypesIndex (const SiteGenotypesIndex& other) : ancestor_index(other.ancestor_index){
//        std::cout << "copy constructor" << std::endl;
//    };
//    SiteGenotypesIndex& operator= (const SiteGenotypesIndex& other) {
//
//        std::cout << "copy assignment" << std::endl;
//        ancestor_index = other.ancestor_index;
//    }

    SiteGenotypesIndex (const SiteGenotypesIndex& other) = default;
    SiteGenotypesIndex& operator= (const SiteGenotypesIndex& other)= default;

//    SiteGenotypesIndex (const SiteGenotypesIndex& other) = delete;
//    SiteGenotypesIndex& operator= (const SiteGenotypesIndex& other) = delete;



    explicit SiteGenotypesIndex(int descendant_count0);

//    ~SiteGenotypesIndex() {}
////        std::cout << "SiteGenotypeIndex destroctor" << std::endl;
//        ancestor_index = 0;
//        descendant_index.clear();
//
//    };

    void SetAncestorIndex(uint32_t index);

    void SetDescendantIndex(int des, uint32_t index);

    uint32_t GetDescendantCount();

    uint32_t GetDescendantIndex(int des);

    uint32_t GetAncestorIndex();

    const std::vector<uint32_t> &GetDescendantIndex();

    void SortIndex();


private:


    uint32_t ancestor_index;
    std::vector<uint32_t> descendant_index ;
};


#endif //_SITE_GENOTYPES_INDEX_H_
