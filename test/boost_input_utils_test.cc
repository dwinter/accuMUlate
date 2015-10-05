//
// Created by steven on 9/30/15.
//
#include <gtest/gtest.h>
#include "unittest_utils.h"
#include "boost_input_utils.h"


namespace po = boost::program_options;
class BoostInputUtilsTest : public ::testing::Test {


public:

protected:

    SamHeader header;
    virtual void SetUp() {
        std::cerr.setstate(std::ios_base::failbit) ; //Supress warning

        header = SamHeader(
                "@RG\tID:M0_CGATGT_L001\tPL:ILLUMINA\tPU:D0V1CACXX.1.TtM0\tLB:LIB-TtM0\tPI:363\tSM:TtM0\tCN:DNASU\n"
                "@RG\tID:M0_CGATGT_L005\tPL:ILLUMINA\tPU:D0V1CACXX.5.TtM0\tLB:LIB-TtM0\tPI:363\tSM:TtM0\tCN:DNASU\n"
                "@RG\tID:M19_AACGTGAT_L007\tPL:ILLUMINA\tPU:C33G8ACXX.7.TtM19\tLB:LIB-TtM19\tPI:350\tSM:TtM19\tCN:IU\n"
                "@RG\tID:M20_AAACATCG_L007\tPL:ILLUMINA\tPU:C33G8ACXX.7.TtM20\tLB:LIB-TtM20\tPI:350\tSM:TtM20\tCN:IU\n"
                "@RG\tID:M25_GTGTTCTA_L007\tPL:ILLUMINA\tPU:C33G8ACXX.7.TtM25\tLB:LIB-TtM25\tPI:350\tSM:TtM25\tCN:IU\n"
                "@RG\tID:M28_TGACCA_L001\tPL:ILLUMINA\tPU:D0V1CACXX.1.TtM28\tLB:LIB-TtM28\tPI:363\tSM:TtM28\tCN:DNASU\n"
                "@RG\tID:M28_TGACCA_L005\tPL:ILLUMINA\tPU:D0V1CACXX.5.TtM28\tLB:LIB-TtM28\tPI:363\tSM:TtM28\tCN:DNASU\n"
                "@RG\tID:M29_CAGATCTG_L007\tPL:ILLUMINA\tPU:C33G8ACXX.7.TtM29\tLB:LIB-TtM29\tPI:350\tSM:TtM29\tCN:IU\n"
                "@RG\tID:M40_ACAGTG_L001\tPL:ILLUMINA\tPU:D0V1CACXX.1.TtM40\tLB:LIB-TtM40\tPI:363\tSM:TtM40\tCN:DNASU\n"
                "@RG\tID:M40_ACAGTG_L005\tPL:ILLUMINA\tPU:D0V1CACXX.5.TtM40\tLB:LIB-TtM40\tPI:363\tSM:TtM40\tCN:DNASU\n"
                "@RG\tID:M44_GCCAAT_L001\tPL:ILLUMINA\tPU:D0V1CACXX.1.TtM44\tLB:LIB-TtM44\tPI:363\tSM:TtM44\tCN:DNASU\n"
                "@RG\tID:M44_GCCAAT_L005\tPL:ILLUMINA\tPU:D0V1CACXX.5.TtM44\tLB:LIB-TtM44\tPI:363\tSM:TtM44\tCN:DNASU\n"
                "@RG\tID:M47_ATTGAGGA_L007\tPL:ILLUMINA\tPU:C33G8ACXX.7.TtM47\tLB:LIB-TtM47\tPI:350\tSM:TtM47\tCN:IU\n"
                "@RG\tID:M50_CAGATC_L001\tPL:ILLUMINA\tPU:D0V1CACXX.1.TtM50\tLB:LIB-TtM50\tPI:363\tSM:TtM50\tCN:DNASU\n"
                "@RG\tID:M50_CAGATC_L005\tPL:ILLUMINA\tPU:D0V1CACXX.5.TtM50\tLB:LIB-TtM50\tPI:363\tSM:TtM50\tCN:DNASU\n"
                "@RG\tID:M51_GGAGAACA_L007\tPL:ILLUMINA\tPU:C33G8ACXX.7.TtM51\tLB:LIB-TtM51\tPI:350\tSM:TtM51\tCN:IU\n"
                "@RG\tID:M531_CTTGTA_L001\tPL:ILLUMINA\tPU:D0V1CACXX.1.TtM531\tLB:LIB-TtM531\tPI:363\tSM:TtM531\tCN:DNASU\n"
                "@RG\tID:M531_CTTGTA_L005\tPL:ILLUMINA\tPU:D0V1CACXX.5.TtM531\tLB:LIB-TtM531\tPI:363\tSM:TtM531\tCN:DNASU\n"
        );

    };

    virtual void TearDown(){
        std::cerr.clear() ; // Add std::cerr back
    }
};



TEST_F(BoostInputUtilsTest, TestSampleMap) {

    std::vector<string> keepers =  {"TtM19", "TtM20", "TtM25", "TtM28", "TtM29",
                                    "TtM40", "TtM44", "TtM47", "TtM50", "TtM51",
                                    "TtM531"};
    string anc_tag = "TtM0";

    po::variables_map vm;
    vm.insert(std::make_pair("sample-name", po::variable_value(keepers, false)));
    vm.insert(std::make_pair("ancestor", po::variable_value(anc_tag, false)));
    po::notify(vm);

    SampleMap samples = BoostUtils::ParseSamples(vm, header);

    SampleMap expected_map;
    expected_map["M0_CGATGT_L001"] = 0;
    expected_map["M0_CGATGT_L005"] = 0;
    expected_map["M19_AACGTGAT_L007"] = 1;
    expected_map["M20_AAACATCG_L007"] = 2;
    expected_map["M25_GTGTTCTA_L007"] = 3;
    expected_map["M28_TGACCA_L001"] = 4;
    expected_map["M28_TGACCA_L005"] = 4;
    expected_map["M29_CAGATCTG_L007"] = 5;
    expected_map["M40_ACAGTG_L001"] = 6;
    expected_map["M40_ACAGTG_L005"] = 6;
    expected_map["M44_GCCAAT_L001"] = 7;
    expected_map["M44_GCCAAT_L005"] = 7;
    expected_map["M47_ATTGAGGA_L007"] = 8;
    expected_map["M50_CAGATC_L001"] = 9;
    expected_map["M50_CAGATC_L005"] = 9;
    expected_map["M51_GGAGAACA_L007"] = 10;
    expected_map["M531_CTTGTA_L001"] = 11;
    expected_map["M531_CTTGTA_L005"] = 11;

    ASSERT_EQ(expected_map.size(), samples.size());
    for (auto item : expected_map) {
        ASSERT_EQ(item.second, samples[item.first]);
    }

}


TEST_F(BoostInputUtilsTest, TestSampleMapNoDup) {

    std::vector<string> keepers =  {"TtM19", "TtM20", "TtM25",
                                    "TtM29",
                                    "TtM47", "TtM51"};
    string anc_tag = "TtM0";

    po::variables_map vm;
    vm.insert(std::make_pair("sample-name", po::variable_value(keepers, false)));
    vm.insert(std::make_pair("ancestor", po::variable_value(anc_tag, false)));
    po::notify(vm);

    SampleMap samples = BoostUtils::ParseSamples(vm, header);

    SampleMap expected_map;
    expected_map["M0_CGATGT_L001"] = 0;
    expected_map["M0_CGATGT_L005"] = 0;
    expected_map["M19_AACGTGAT_L007"] = 1;
    expected_map["M20_AAACATCG_L007"] = 2;
    expected_map["M25_GTGTTCTA_L007"] = 3;
    expected_map["M28_TGACCA_L001"] = std::numeric_limits<uint32_t>::max();
    expected_map["M28_TGACCA_L005"] = std::numeric_limits<uint32_t>::max();
    expected_map["M29_CAGATCTG_L007"] = 4;
    expected_map["M40_ACAGTG_L001"] = std::numeric_limits<uint32_t>::max();
    expected_map["M40_ACAGTG_L005"] = std::numeric_limits<uint32_t>::max();
    expected_map["M44_GCCAAT_L001"] = std::numeric_limits<uint32_t>::max();
    expected_map["M44_GCCAAT_L005"] = std::numeric_limits<uint32_t>::max();
    expected_map["M47_ATTGAGGA_L007"] = 5;
    expected_map["M50_CAGATC_L001"] = std::numeric_limits<uint32_t>::max();
    expected_map["M50_CAGATC_L005"] = std::numeric_limits<uint32_t>::max();
    expected_map["M51_GGAGAACA_L007"] = 6;
    expected_map["M531_CTTGTA_L001"] = std::numeric_limits<uint32_t>::max();
    expected_map["M531_CTTGTA_L005"] = std::numeric_limits<uint32_t>::max();

    ASSERT_EQ(expected_map.size(), samples.size());
    for (auto item : expected_map) {
        ASSERT_EQ(item.second, samples[item.first]);
    }

}


TEST_F(BoostInputUtilsTest, TestSampleMapOnlyDup) {

    std::vector<string> keepers =  {"TtM28",
                                    "TtM40", "TtM44",  "TtM50",
                                    "TtM531"};

    string anc_tag = "TtM0";

    po::variables_map vm;
    vm.insert(std::make_pair("sample-name", po::variable_value(keepers, false)));
    vm.insert(std::make_pair("ancestor", po::variable_value(anc_tag, false)));
    po::notify(vm);

    SampleMap samples = BoostUtils::ParseSamples(vm, header);

    SampleMap expected_map;
    expected_map["M0_CGATGT_L001"] = 0;
    expected_map["M0_CGATGT_L005"] = 0;
    expected_map["M19_AACGTGAT_L007"] = std::numeric_limits<uint32_t>::max();
    expected_map["M20_AAACATCG_L007"] = std::numeric_limits<uint32_t>::max();
    expected_map["M25_GTGTTCTA_L007"] = std::numeric_limits<uint32_t>::max();
    expected_map["M28_TGACCA_L001"] = 1;
    expected_map["M28_TGACCA_L005"] = 1;
    expected_map["M29_CAGATCTG_L007"] = std::numeric_limits<uint32_t>::max();
    expected_map["M40_ACAGTG_L001"] = 2;
    expected_map["M40_ACAGTG_L005"] = 2;
    expected_map["M44_GCCAAT_L001"] = 3;
    expected_map["M44_GCCAAT_L005"] = 3;
    expected_map["M47_ATTGAGGA_L007"] = std::numeric_limits<uint32_t>::max();
    expected_map["M50_CAGATC_L001"] = 4;
    expected_map["M50_CAGATC_L005"] = 4;
    expected_map["M51_GGAGAACA_L007"] = std::numeric_limits<uint32_t>::max();
    expected_map["M531_CTTGTA_L001"] = 5;
    expected_map["M531_CTTGTA_L005"] = 5;

    ASSERT_EQ(expected_map.size(), samples.size());
    for (auto item : expected_map) {
        ASSERT_EQ(item.second, samples[item.first]);
    }

}


TEST_F(BoostInputUtilsTest, TestSampleMapRandom1) {

    std::vector<string> keepers =  {"TtM19", "TtM20", "TtM25", "TtM28", "TtM29",
//                                    "TtM40", "TtM44", "TtM47", "TtM50", "TtM51",
                                    "TtM531"};

    string anc_tag = "TtM0";

    po::variables_map vm;
    vm.insert(std::make_pair("sample-name", po::variable_value(keepers, false)));
    vm.insert(std::make_pair("ancestor", po::variable_value(anc_tag, false)));
    po::notify(vm);

    SampleMap samples = BoostUtils::ParseSamples(vm, header);

    SampleMap expected_map;
    expected_map["M0_CGATGT_L001"] = 0;
    expected_map["M0_CGATGT_L005"] = 0;
    expected_map["M19_AACGTGAT_L007"] = 1;
    expected_map["M20_AAACATCG_L007"] = 2;
    expected_map["M25_GTGTTCTA_L007"] = 3;
    expected_map["M28_TGACCA_L001"] = 4;
    expected_map["M28_TGACCA_L005"] = 4;
    expected_map["M29_CAGATCTG_L007"] = 5;
    expected_map["M40_ACAGTG_L001"] = std::numeric_limits<uint32_t>::max();
    expected_map["M40_ACAGTG_L005"] = std::numeric_limits<uint32_t>::max();
    expected_map["M44_GCCAAT_L001"] = std::numeric_limits<uint32_t>::max();
    expected_map["M44_GCCAAT_L005"] = std::numeric_limits<uint32_t>::max();
    expected_map["M47_ATTGAGGA_L007"] = std::numeric_limits<uint32_t>::max();
    expected_map["M50_CAGATC_L001"] = std::numeric_limits<uint32_t>::max();
    expected_map["M50_CAGATC_L005"] = std::numeric_limits<uint32_t>::max();
    expected_map["M51_GGAGAACA_L007"] = std::numeric_limits<uint32_t>::max();
    expected_map["M531_CTTGTA_L001"] = 6;
    expected_map["M531_CTTGTA_L005"] = 6;

    ASSERT_EQ(expected_map.size(), samples.size());
    for (auto item : expected_map) {
        ASSERT_EQ(item.second, samples[item.first]);
    }

}



TEST_F(BoostInputUtilsTest, TestSampleMapRandom2) {

    std::vector<string> keepers =  {
//            "TtM19", "TtM20", "TtM25", "TtM28", "TtM29",
                                    "TtM40", "TtM44", "TtM47", "TtM50", "TtM51"
//                                    "TtM531"
    };

    string anc_tag = "TtM0";

    po::variables_map vm;
    vm.insert(std::make_pair("sample-name", po::variable_value(keepers, false)));
    vm.insert(std::make_pair("ancestor", po::variable_value(anc_tag, false)));
    po::notify(vm);

    SampleMap samples = BoostUtils::ParseSamples(vm, header);

    SampleMap expected_map;
    expected_map["M0_CGATGT_L001"] = 0;
    expected_map["M0_CGATGT_L005"] = 0;
    expected_map["M19_AACGTGAT_L007"] = std::numeric_limits<uint32_t>::max();
    expected_map["M20_AAACATCG_L007"] = std::numeric_limits<uint32_t>::max();
    expected_map["M25_GTGTTCTA_L007"] = std::numeric_limits<uint32_t>::max();
    expected_map["M28_TGACCA_L001"] = std::numeric_limits<uint32_t>::max();
    expected_map["M28_TGACCA_L005"] = std::numeric_limits<uint32_t>::max();
    expected_map["M29_CAGATCTG_L007"] = std::numeric_limits<uint32_t>::max();
    expected_map["M40_ACAGTG_L001"] = 1;
    expected_map["M40_ACAGTG_L005"] = 1;
    expected_map["M44_GCCAAT_L001"] = 2;
    expected_map["M44_GCCAAT_L005"] = 2;
    expected_map["M47_ATTGAGGA_L007"] = 3;
    expected_map["M50_CAGATC_L001"] = 4;
    expected_map["M50_CAGATC_L005"] = 4;
    expected_map["M51_GGAGAACA_L007"] = 5;
    expected_map["M531_CTTGTA_L001"] = std::numeric_limits<uint32_t>::max();
    expected_map["M531_CTTGTA_L005"] = std::numeric_limits<uint32_t>::max();

    ASSERT_EQ(expected_map.size(), samples.size());
    for (auto item : expected_map) {
        ASSERT_EQ(item.second, samples[item.first]);
    }

}


TEST_F(BoostInputUtilsTest, TestAncestorNotExist) {

    std::vector<string> keepers =  {"TtM19"};
    string anc_tag = "TtM";

    po::variables_map vm;
    vm.insert(std::make_pair("sample-name", po::variable_value(keepers, false)));
    vm.insert(std::make_pair("ancestor", po::variable_value(anc_tag, false)));
    po::notify(vm);

    //    SampleMap samples = BoostUtils::ParseSamples(vm, header);
    ASSERT_EXIT(BoostUtils::ParseSamples(vm, header), ::testing::ExitedWithCode(5), "");

}




TEST_F(BoostInputUtilsTest, TestSampleNotExist) {

    std::vector<string> keepers =  {
            "TtM19", "Not_exist"
    };

    string anc_tag = "TtM0";

    po::variables_map vm;
    vm.insert(std::make_pair("sample-name", po::variable_value(keepers, false)));
    vm.insert(std::make_pair("ancestor", po::variable_value(anc_tag, false)));
    po::notify(vm);

//    SampleMap samples = BoostUtils::ParseSamples(vm, header);
    ASSERT_EXIT(BoostUtils::ParseSamples(vm, header), ::testing::ExitedWithCode(6), "");

}


TEST_F(BoostInputUtilsTest, TestCreateModelParams){
    //Simple function, BUT assuming fix order
    //ModelParams at the time of this test is created
//    struct ModelParams{
//        double theta;               //
//        std::vector<double> nuc_freq;    //ACGT
//        double mutation_rate;       //
//        double error_prob;          // Sequencing error-rate
//        double phi_haploid;         // Overdispersion for haploid sequencing
//        double phi_diploid;         // Overdispersion for diploid sequencing
//        int ploidy_ancestor;         //
//        int ploidy_descendant;       //
//    };

    double theta = 0.01;               //
    nfreqs nuc_freq = {{0.25, 0.1, 0.3, 0.35}};    //ACGT
    double mutation_rate = 0.02;       //
    double error_prob = 0.03;          // Sequencing error-rate
    double phi_haploid = 0.04;         // Overdispersion for haploid sequencing
    double phi_diploid = 0.05;         // Overdispersion for diploid sequencing
    int ploidy_ancestor = 1;         //
    int ploidy_descendant = 2;


    po::variables_map vm;
//    vm.insert(std::make_pair("sample-name", po::variable_value(keepers, false)));
//    vm.insert(std::make_pair("ancestor", po::variable_value(anc_tag, false)));
    vm.insert(std::make_pair("theta", po::variable_value(theta, false)));
    vm.insert(std::make_pair("nfreqs", po::variable_value(nuc_freq, false)));
    vm.insert(std::make_pair("mu", po::variable_value(mutation_rate, false)));
    vm.insert(std::make_pair("seq-error", po::variable_value(error_prob, false)));
    vm.insert(std::make_pair("phi-haploid", po::variable_value(phi_haploid, false)));
    vm.insert(std::make_pair("phi-diploid", po::variable_value(phi_diploid, false)));
    vm.insert(std::make_pair("ploidy-ancestor", po::variable_value(ploidy_ancestor, false)));
    vm.insert(std::make_pair("ploidy-descendant", po::variable_value(ploidy_descendant, false)));
    po::notify(vm);
    
    ModelParams model = BoostUtils::CreateModelParams(vm);
    ASSERT_EQ(theta, model.theta);
    ASSERT_EQ(error_prob, model.error_prob);
    ASSERT_EQ(mutation_rate, model.mutation_rate);
    ASSERT_EQ(phi_haploid, model.phi_haploid);
    ASSERT_EQ(phi_diploid, model.phi_diploid);
    ASSERT_EQ(ploidy_ancestor, model.ploidy_ancestor);
    ASSERT_EQ(ploidy_descendant, model.ploidy_descendant);
    for (int i = 0; i < nuc_freq.freqs.size(); ++i) {
        ASSERT_EQ(nuc_freq.freqs[i], model.nuc_freq[i]);
    }
}




// print boost vm
//    for (const auto& it : vm) {
//        std::cout << it.first.c_str() << " ";
//        auto& value = it.second.value();
//        if (auto v = boost::any_cast<int>(&value))
//            std::cout << *v;
//        else if (auto v = boost::any_cast<std::string>(&value))
//            std::cout << *v;
//        else
//            std::cout << "error: " <<  "\t" << &value << "\n";
//    }