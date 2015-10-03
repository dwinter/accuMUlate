//
// Created by steven on 9/30/15.
//
#include <gtest/gtest.h>
#include "variant_visitor.h"


class ParserTest : public ::testing::Test {

public:

protected:

    virtual void SetUp() {

    };
};


TEST_F(ParserTest, TestParserConstructor) {
// Constructors
//    VariantVisitor(const RefVector& bam_references,
//                   LocalBamToolsUtils::Fasta& idx_ref,
//                   ostream *out_stream,
//                   SampleMap& samples,
//                   const ModelParams& p,
//                   BamAlignment& ali,
//                   int qual_cut, int mapping_cut, double prob_cut);
//    ReadDataVisitor(const RefVector& bam_references,
//                    LocalBamToolsUtils::Fasta& idx_ref,
//                    SampleMap& samples,
//                    const ModelParams& p,
//                    BamAlignment& ali,
//                    int qual_cut,
//                    int mapping_cut);

    streambuf * buf;

    RefVector references;
    LocalBamToolsUtils::Fasta reference_genome; // BamTools::Fasta
    ostream result_stream(buf);
    SampleMap samples;
    ModelParams params;
    BamAlignment ali;

    int exp_qual = 10;
    int exp_map_qual = 20;
    double prob = 0.1;
//    VariantVisitor v(
//            references,
//            reference_genome,
//            &result_stream,
//            samples,
//            params,
//            ali,
//            10,20,0.1
//    );

}


TEST_F(ParserTest, Testbase_index_lookup) {


    ASSERT_EQ(base_index_lookup['a'], 0);
    ASSERT_EQ(base_index_lookup['A'], 0);
    ASSERT_EQ(base_index_lookup['c'], 1);
    ASSERT_EQ(base_index_lookup['C'], 1);
    ASSERT_EQ(base_index_lookup['g'], 2);
    ASSERT_EQ(base_index_lookup['G'], 2);
    ASSERT_EQ(base_index_lookup['t'], 3);
    ASSERT_EQ(base_index_lookup['T'], 3);

}