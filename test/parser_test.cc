//
// Created by steven on 9/30/15.
//
#include <gtest/gtest.h>
#include "variant_visitor.h"

#include "boost/program_options.hpp"
#include "api/BamReader.h"
#include "boost_input_utils.h"

#include "src/io_data/local_bamtools/bamtools_pileup_engine.h"
#include "src/io_data/local_bamtools/bamtools_fasta.h"



class ParserTest : public ::testing::Test {

public:
    ParserTest() {}

        BamAlignment ali;
        BamReader experiment;
        SamHeader header;
        RefVector references;
        LocalBamToolsUtils::Fasta reference_genome; // BamTools::Fasta
        SampleMap samples;
        vector<string> sample_names = {"A0", "D1", "D2", "D3", "D4"};
        ModelParams params = { 
            0.0001, 
            {0.38, 0.12, 0.12, 0.38}, 
            1e-8,
            0.01,
            0.001,
            0.005,
            2,2
        };
 



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
    ostream result_stream(buf);

    VariantVisitor v(
            references,
            reference_genome,
            &result_stream,
            samples,
            sample_names,
            params,
            ali,
            10,20,0.1
    );

}

TEST_F(ParserTest, parse_bed) {
    BedFile bed ("data/test.bed");
    BedInterval region;
    int i = 0;
    while(bed.get_interval(region) == 0){        
        i++;        
    }
    ASSERT_EQ(i, 3);
    ASSERT_EQ(region.chr, "double_mutation");
    ASSERT_EQ(region.start, 201);
    ASSERT_EQ(region.end, 203);
}

TEST_F(ParserTest, restrict_to_interval){

    //Do we _only_ call mutations from within a given interval
    //when that interval is specified. Bugs #11 and #28 lead to the progarm
    //calling form surrounding intervals and only the final interval in a file.
    //Rather than capture ouptur or dig into the visitor, we generate an output
    //file and and parse that (as BED) to confirm the correct behaviour


    // Unfortunately, doing this without actually parsing a command line seems
    // to take a lot of copy-pasting...
    
    namespace po = boost::program_options;  

    po::options_description cmd("Command line options");
    cmd.add_options()
        ("help,h", "Print a help message")
        ("bam,b", po::value<string>()->required(), "Path to BAM file")
        ("bam-index,x", po::value<string>()->default_value(""), "Path to BAM index, (defalult is <bam_path>.bai")
        ("reference,r", po::value<string>()->required(),  "Path to reference genome")
        ("ancestor,a", po::value<string>()->required(), "Ancestor RG sample ID")
        ("sample-name,s", po::value<vector <string> >()->required(), "Sample tags to include")
        ("intervals,i", po::value<string>(), "Path to bed file")
        ("config,c", po::value<string>(), "Path to config file")
        ("theta", po::value<double>()->required(), "theta")
//        ("nfreqs", po::value< nfreqs >()->multitoken(), "Nucleotide frequencies")
        ("mu", po::value<double>()->required(), "Experiment-long mutation rate")
        ("seq-error", po::value<double>()->required(), "Probability of sequencing error")
        ("ploidy-ancestor", po::value<int>()->default_value(2), "Polidy of ancestor (1 or 2)")
        ("ploidy-descendant", po::value<int>()->default_value(2), "Ploidy of descendant (1 or 2)")
        ("phi-haploid",     po::value<double>(), "Over-dispersion for haploid sequencing")
        ("phi-diploid",     po::value<double>(), "Over-dispersion for diploid sequencing");


    boost::program_options::variables_map vm;
    ifstream config_stream("data/test_params.ini");
    po::store(po::parse_config_file(config_stream, cmd, false), vm);
    vm.notify();
    BoostUtils::ExtractInputVariables(vm, experiment, references, header, reference_genome);  

    SampleMap interval_samples = BoostUtils::ParseSamples(vm, header);
    LocalBamToolsUtils::PileupEngine pileup;


    //convoluted, but wriiten to match/test procedure in main.cc 
    streambuf * buf;
    ofstream of;
    of.open("data/test.out");
    buf = of.rdbuf();
    ostream result_stream(buf);

    VariantVisitor *v = new VariantVisitor(
            references,
            reference_genome,
            &result_stream,
            interval_samples,
            sample_names,
            params,
            ali,
            10,20,0.0

    );
    pileup.AddVisitor(v);
    BedFile bed (vm["intervals"].as<string>());
    BedInterval region;
    int ref_id = -1;
    while(bed.get_interval(region) == 0){        
        ref_id = experiment.GetReferenceID(region.chr);
        experiment.SetRegion(ref_id, region.start, ref_id, region.end); 
        while( experiment.GetNextAlignment(ali) ){
                pileup.AddAlignment(ali);
        }
    }

    BedFile result ("data/test.out");
    int i = 0;
    while(bed.get_interval(region) == 0){        
        i++;        
    }

    ASSERT_EQ(i, 4);
    ASSERT_EQ(region.chr, "double_mutation");
    ASSERT_EQ(region.start, 202);
    ASSERT_EQ(region.end, 203);
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
