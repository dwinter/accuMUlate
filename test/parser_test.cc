//
// Created by steven on 9/30/15.
//
#include <gtest/gtest.h>
#include "parsers.h"
#include "unittest_utils.h"



const int BASE_COUNT = 4;

class ParserTest : public ::testing::Test {

public:

protected:

    virtual void SetUp() {

    };
};


//TEST_F(ParserTest, TestParser) {
//
//    FastaReference reference_g("test/test.fai");
//    int chr_idx;
//    reference_g.get_ref_id("scf_8254727", chr_idx);
//    //    cout << "chr index = " << chr_idx << " (should be 9)" << endl;
//    ASSERT_EQ(9, chr_idx);
//
//    BedFile bed("test/test.bed");
//    BedInterval current_line;
//    while (bed.get_interval(current_line) == 0) {
//        cout << current_line.chr << endl;
//    }
//}

