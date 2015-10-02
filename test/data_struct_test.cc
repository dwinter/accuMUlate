//
// Created by steven on 4/8/15.
//
#include <gtest/gtest.h>
#include "data_struct.h"
#include "unittest_utils.h"

const int BASE_COUNT = 4;


class DataStructTest : public ::testing::Test {

public:

protected:

    virtual void SetUp() {}
};


TEST_F(DataStructTest, TestReadData) {

    ReadData rdata = {0};
    ASSERT_EQ(0, rdata.key);
    for (int i = 0; i < BASE_COUNT; ++i) {
        ASSERT_EQ(0, rdata.reads[i]);
    }

    rdata = {10};
    ASSERT_EQ(10, rdata.key);
    for (int i = 1; i < BASE_COUNT; ++i) {
        ASSERT_EQ(0, rdata.reads[i]);
    }
    ASSERT_EQ(10, rdata.reads[0]);

    rdata = {65535};
    ASSERT_EQ(65535, rdata.key);
    for (int i = 1; i < BASE_COUNT; ++i) {
        ASSERT_EQ(0, rdata.reads[i]);
    }
    ASSERT_EQ(65535, rdata.reads[0]);

    rdata = {65536};
    ASSERT_EQ(65536, rdata.key);
    for (int i : {0, 2, 3}) {
        ASSERT_EQ(0, rdata.reads[i]);
    }
    ASSERT_EQ(1, rdata.reads[1]);

    rdata = ReadData(0);
    rdata.reads[2] = 1;
    ASSERT_EQ(4294967296, rdata.key);

    rdata = ReadData{281474976710656};
    for (int i : {0, 1, 2}) {
        ASSERT_EQ(0, rdata.reads[i]);
    }
    ASSERT_EQ(1, rdata.reads[3]);

//    281474976710656 + 4294967296 + 65536 + 1
//    [1] 281479271743489
    rdata = {281479271743489};
    for (int i : {0, 1, 2, 3}) {
        ASSERT_EQ(1, rdata.reads[i]);
    }

//    281474976710656*5 + 4294967296*4 + 65536*3 + 1*2
//    [1] 1407392063619074
    uint16_t temp_r[4] = {2, 3, 4, 5}; // array_var
	rdata = ReadData {temp_r};
    for (int i : {0, 1, 2, 3}) {
        ASSERT_EQ(i+2, rdata.reads[i]);
    }
    ASSERT_EQ(1407392063619074, rdata.key);

	ReadData r4 ({2, 3, 4, 5} ); // initializer_list reads
    ASSERT_EQ(1407392063619074, r4.key);
	ReadData r5 {{2, 3, 4, 5} }; // initializer_list reads
    ASSERT_EQ(1407392063619074, r5.key);


}