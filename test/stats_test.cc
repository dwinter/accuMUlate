#include <gtest/gtest.h>
#include <src/stats/stats.h>
#include "unittest_utils.h"


class StatsTest : public ::testing::Test {
    public:

    protected:
        
};
    

TEST_F (StatsTest, fisher_few_reads) {
    //Really the exact test
    ASSERT_NEAR(fisher_exact_test(6,12,12,5), 0.0437102, 1e-7);
}

TEST_F (StatsTest, fisher_many_reads){
    //A gtest because there are too many reads to do the exact one.
    //known p-value is from a simulation, so error margin is lowered.
    ASSERT_NEAR(fisher_exact_test(32,18,241,20), 9.12e-07, 1e-7);
}

TEST_F(StatsTest, ad_test){
    //tested against scipy k-sample test
    std::vector<int> a = {40, 31, 35, 40, 40, 32, 33};           
    std::vector<int> b = {21, 31, 33, 34, 34, 40, 42, 20} ;
    ASSERT_NEAR(ad_two_sample_test(a,b), -0.52579911592960638, ERROR_THRESHOLD);
}

TEST_F(StatsTest, phred){
    ASSERT_NEAR(phred(0.05), 13.0103, 1e-4);
}


