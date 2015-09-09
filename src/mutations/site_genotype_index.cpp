/*
 * site_genotypes_index.cc
 *
 *  Created on: 4/5/15
 *      Author: Steven Wu
 */

//
// Created by steven on 4/5/15.
//

#include "site_genotype_index.h"

/*
 * sequence_prob.cc
 *
 *  Created on: Nov 7, 2014
 *      Author: Steven Wu
 */

//SiteGenotypesIndex::~SiteGenotypesIndex() {
//}


SiteGenotypesIndex::SiteGenotypesIndex(int descendant_count0) {

    descendant_index.assign(descendant_count0, -1);
}

void SiteGenotypesIndex::SetAncestorIndex(uint32_t index) {
    ancestor_index = index;
}

void SiteGenotypesIndex::SetDescendantIndex(int des, uint32_t index) {
    descendant_index[des] = index;
}


uint32_t SiteGenotypesIndex::GetDescendantCount() {
    return descendant_index.size();
}


uint32_t SiteGenotypesIndex::GetDescendantIndex(int des) {
    return descendant_index[des];

}

const std::vector<uint32_t> &SiteGenotypesIndex::GetDescendantIndex() {
    return descendant_index;
}


uint32_t SiteGenotypesIndex::GetAncestorIndex() {
    return ancestor_index;
}

void SiteGenotypesIndex::SortIndex() {
//	descendant_index
//	std::algorithm::sort
    std::sort(descendant_index.begin(), descendant_index.end());

}
