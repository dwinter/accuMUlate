/*
 * boost_input_utils.cc
 *
 *  Created on: 12/12/14
 *      Author: Steven Wu
 *  Modified DJW Aug-15
 */


#include "parsers.h"
#include "boost_input_utils.h"

void validate(boost::any& v, const vector<string>& values, nfreqs* target_type, int){
    nfreqs result;
    vector<double> pi;
    for(vector<string>::const_iterator it = values.begin(); 
        it != values.end(); ++it){
        stringstream ss(*it);
        copy(istream_iterator<double>(ss), istream_iterator<double>(), back_inserter(pi));
    }
    if(pi.size() != 4){
        throw boost::program_options::invalid_option_value("Must specify 4 (and only 4) nucleotide frequencies");        
    }
    if( accumulate(pi.begin(), pi.end(), 0.0) != 1.0){
        throw boost::program_options::invalid_option_value("Nucleotide frequencies don't sum to 1.0");          }
    result.freqs = pi;
    v= result;
}

namespace BoostUtils {
    using namespace std;
    namespace po = boost::program_options;

    static bool file_exists(const std::string &name) {
        struct stat buffer;
        return (stat(name.c_str(), &buffer) == 0);
    }
    
    SampleMap ParseSamples(boost::program_options::variables_map &vm, BamTools::SamHeader &header){
        //OK. We want to store data form each sample in a unique position in 
        //a ReadDataVector. The caller also assumes that the ancestral sequence
        //is in the 0th poistion of the ReadDataVector.
        //The sample info in the BAM file is all in the ReadGroup headers, and
        //we want to be able to exclude some samples.So, we ask users to provide
        //a list of included samples.
        //To acheive all this we start parsing through every read group in the
        //BAM. If it's one to include, we add it the index map and remove it
        //from the list of samples to include. If it's not in our include list
        //we warn the user we are skipping some of the data. 
        SampleMap name_map;
        vector<string> keepers =  vm["sample-name"].as< vector<string> >();
        for (auto key : keepers) {
            cout << "sample-name: " << key << endl;
        }
        string anc_tag = vm["ancestor"].as<string>();
        cout << "A: " << anc_tag << endl;
        uint16_t sindex = 1;
        bool ancestor_in_BAM = false;
        for(auto it = header.ReadGroups.Begin(); it!= header.ReadGroups.End(); it++){
            if(it->HasSample()){

                if (find(keepers.begin(), keepers.end(), it->Sample) == keepers.end()){
                    if(it->Sample == anc_tag ){
                          name_map[it->Sample] = 0;
                          ancestor_in_BAM = true;
                    }
                    else {
                        name_map[it->Sample] = numeric_limits<uint32_t>::max()  ;
                        cerr << "Warning: excluding data from '" << it->Sample <<
                             "' which is included in the BAM file but not the list of included samples" << endl;
                    }
                }
                else {
                    auto s  = name_map.find(it->Sample);
                    if( s == name_map.end()){//samples can have multiple readgroups...
                        name_map[it->Sample] = sindex;
                        sindex += 1;
//                        keepers.erase(find(keepers.begin(),keepers.end(),it->Sample));
                        //NOTE: with erase, effectively remove all double+ samples
                        //Fixed version == #sample-name M28 40 44 50 531
                    }
                }

            }
        }
//Original:
//        12	0	12
//        TtM531	4294967295
//        TtM51	10
//        TtM25	3
//        TtM47	8
//        TtM28	4294967295
//        TtM19	1
//        TtM0	0
//        TtM29	5
//        TtM40	4294967295
//        TtM20	2
//        TtM44	4294967295
//        TtM50	4294967295
//        M0_CGATGT_L001 0
//        M0_CGATGT_L005 0
//        M19_AACGTGAT_L007 1
//        M20_AAACATCG_L007 2
//        M25_GTGTTCTA_L007 3
//        M28_TGACCA_L001 4294967295
//        M28_TGACCA_L005 4294967295
//        M29_CAGATCTG_L007 5
//        M40_ACAGTG_L001 4294967295
//        M40_ACAGTG_L005 4294967295
//        M44_GCCAAT_L001 4294967295
//        M44_GCCAAT_L005 4294967295
//        M47_ATTGAGGA_L007 8
//        M50_CAGATC_L001 4294967295
//        M50_CAGATC_L005 4294967295
//        M51_GGAGAACA_L007 10
//        M531_CTTGTA_L001 4294967295
//        M531_CTTGTA_L005 4294967295
// update version without erase
//        TtM531	4294967295
//        TtM51	6
//        TtM25	3
//        TtM47	5
//        TtM28	4294967295
//        TtM19	1
//        TtM0	0
//        TtM29	4
//        TtM40	4294967295
//        TtM20	2
//        TtM44	4294967295
//        TtM50	4294967295
//        M0_CGATGT_L001 0
//        M0_CGATGT_L005 0
//        M19_AACGTGAT_L007 1
//        M20_AAACATCG_L007 2
//        M25_GTGTTCTA_L007 3
//        M28_TGACCA_L001 4294967295
//        M28_TGACCA_L005 4294967295
//        M29_CAGATCTG_L007 4
//        M40_ACAGTG_L001 4294967295
//        M40_ACAGTG_L005 4294967295
//        M44_GCCAAT_L001 4294967295
//        M44_GCCAAT_L005 4294967295
//        M47_ATTGAGGA_L007 5
//        M50_CAGATC_L001 4294967295
//        M50_CAGATC_L005 4294967295
//        M51_GGAGAACA_L007 6
//        M531_CTTGTA_L001 4294967295
//        M531_CTTGTA_L005 4294967295

        //Once we've built the index map we can check if we now about every
        //sample in the BAM file and if we have set the ancesoe
        if(!ancestor_in_BAM){
            cerr << "Error: No data for ancestral sample '" << anc_tag <<
                    "' in the specifified BAM file. Check the sample tags match" << endl;
            exit(1);
        }
        cout << name_map.size() << "\t" << keepers.size() << "\t" << sindex << endl;
        for (auto item : name_map) {
            cout << item.first << "\t" << item.second << endl;
        }
        if( (keepers.size()+1) != sindex ){
            cerr << "Sample(s) note persent in BAM file: ";
            for(auto s: keepers){
                cerr << s << " ";
            }
            cerr << endl;
//            exit(1);
        }
        // And now.. go back over the read groups to map RG->sample index
        SampleMap samples;
        for(auto it = header.ReadGroups.Begin(); it!= header.ReadGroups.End(); it++){
            if(it->HasSample()){
                samples[it->ID] = name_map[it->Sample];
                cout << it->ID << " " << name_map[it->Sample] << endl;
            }
        }


//
//
//
////EM version
//        SampleMap sample_map_temp;
//        int total_sample_count = 0;
//        SamReadGroupDictionary dictionary = header.ReadGroups;
//        size_t tag_count = dictionary.Size();
//
//        for (auto it = dictionary.Begin(); it != dictionary.End(); it++) {
//            if (it->HasSample()) {
//                auto s = sample_map_temp.find(it->Sample);
//                if (s == sample_map_temp.end()) { // not in there yet
//                    sample_map_temp[it->Sample] = total_sample_count;
//                    total_sample_count++;
//                }
//            }
//        }
//
////        m_samples = sample_map_temp;
//
//        SampleMap map_tag_sample = std::unordered_map<std::string, uint32_t>(tag_count);
//        for (auto dict = dictionary.Begin(); dict != dictionary.End(); dict++) {
//            map_tag_sample.emplace(dict->ID, sample_map_temp[dict->Sample]);
////        map_tag_sample_two_stage.emplace(dict->ID, dict->Sample);
//        }
//        std::cout << "================= V2 =============" << std::endl;
//        for (auto sample : map_tag_sample) {
//            std::cout << sample.first << "\t" << sample.second << std::endl;
//        }
//
//

        return samples;

    }


    
    void check_args(boost::program_options::variables_map &vm){
        // Is the experimental design one of the ones we can handle?
        int ploidy_ancestor = vm["ploidy-ancestor"].as<int>();
        int polidy_descendant = vm["ploidy-descendant"].as<int>();

        if (ploidy_ancestor > 2 or ploidy_ancestor < 1){
            throw po::invalid_option_value("accuMUlate can't only deal with haploid or diploid ancestral samples"); 
        }
        if (polidy_descendant > 2 or polidy_descendant < 1){
            throw po::invalid_option_value("accuMUlate can't only deal with haploid or descendant samples");        
        }
        if (ploidy_ancestor == 1 and polidy_descendant == 2){
            throw po::invalid_option_value("accuMUlate has no model for a haploid->diploid MA experiemt");
        }
        //Do we have the right over-dispersion params set
        if (ploidy_ancestor == 1 or polidy_descendant == 1){
            if(not vm.count("phi-haploid")){
                throw po::invalid_option_value("Must specify phi-haploid (overdispersion for haploid sequencing)");
            }
        }
        if (ploidy_ancestor == 2 or polidy_descendant == 2){
            if(not vm.count("phi-diploid")){
                throw po::invalid_option_value("Must specify phi-diploid (overdispersion for diploid sequencing)");
            }
        }

    }


    void ParseCommandLineInput(int argc, char **argv, boost::program_options::variables_map &vm) {
        po::options_description cmd("Command line options");
        cmd.add_options()
            ("help,h", "Print a help message")
            ("bam,b", po::value<string>()->required(), "Path to BAM file")
            ("bam-index,x", po::value<string>()->default_value(""), "Path to BAM index, (defalult is <bam_path>.bai")
            ("reference,r", po::value<string>()->required(),  "Path to reference genome")
            ("ancestor,a", po::value<string>()->required(), "Ancestor RG sample ID")
            ("sample-name,s", po::value<vector <string> >()->required(), "Sample tags to include")
            ("qual,q", po::value<int>()->default_value(13), "Base quality cuttoff")
            ("mapping-qual,m", po::value<int>()->default_value(13), "Mapping quality cuttoff")
            ("prob,p", po::value<double>()->default_value(0.1), "Mutaton probability cut-off")
            ("out,o", po::value<string>()->default_value(""), "Out file name (default is std out)")
            ("intervals,i", po::value<string>(), "Path to bed file")
            ("config,c", po::value<string>(), "Path to config file")
            ("theta", po::value<double>()->required(), "theta")
            ("nfreqs", po::value< nfreqs >()->multitoken(), "Nucleotide frequencies")
            ("mu", po::value<double>()->required(), "Experiment-long mutation rate")
            ("seq-error", po::value<double>()->required(), "Probability of sequencing error")
            ("ploidy-ancestor", po::value<int>()->default_value(2), "Polidy of ancestor (1 or 2)")
            ("ploidy-descendant", po::value<int>()->default_value(2), "Ploidy of descendant (1 or 2)")
            ("phi-haploid",     po::value<double>(), "Over-dispersion for haploid sequencing")
            ("phi-diploid",     po::value<double>(), "Over-dispersion for diploid sequencing");
        po::store(po::parse_command_line(argc, argv, cmd), vm);
    


        if (vm.count("help")) {
            cout << cmd << endl;
            exit(-1);
        }

        if (vm.count("config")) {
            ifstream config_stream(vm["config"].as<string>());
            po::store(po::parse_config_file(config_stream, cmd, false), vm);
        }

        vm.notify();
        check_args(vm);

    }

    // Set up everything that has to be refered to by reference
    void ExtractInputVariables(boost::program_options::variables_map &vm,
            BamTools::BamReader &experiment, BamTools::RefVector &references,
            BamTools::SamHeader &header, LocalBamToolsUtils::Fasta &reference_genome){


        string ref_file = vm["reference"].as<string>();
        string bam_path = vm["bam"].as<string>();
        string index_path = vm["bam-index"].as<string>();
        if (index_path == "") {
            index_path = bam_path + ".bai";
        }
        experiment.Open(bam_path);
        experiment.OpenIndex(index_path);
        references = experiment.GetReferenceData();
        header = experiment.GetHeader();

        experiment.OpenIndex(index_path);


        if (!file_exists(ref_file + ".fai")) {
            reference_genome.Open(ref_file);
            reference_genome.CreateIndex(ref_file + ".fai");
        }
        else {
            reference_genome.Open(ref_file, ref_file + ".fai");
        }
    }


    ModelParams CreateModelParams(boost::program_options::variables_map vm) {
        ModelParams params = {
            vm["theta"].as<double>(),
            vm["nfreqs"].as< nfreqs >().freqs,
            vm["mu"].as<double>(),
            vm["seq-error"].as<double>(),
            vm["phi-haploid"].as<double>(),
            vm["phi-diploid"].as<double>(),
            vm["ploidy-ancestor"].as<int>(),
            vm["ploidy-descendant"].as<int>()
        };
        return params;

    }
}
