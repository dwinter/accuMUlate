/*
 * boost_input_utils.cc
 *
 *  Created on: 12/12/14
 *      Author: Steven Wu
 */


#include "boost_input_utils.h"
namespace BoostUtils {
    using namespace std;

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
        string anc_tag = vm["ancestor"].as<string>();
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
                        keepers.erase(find(keepers.begin(),keepers.end(),it->Sample));
                    }
                }
            }
        }
        //Once we've built the index map we can check if we now about every
        //sample in the BAM file and if we have set the ancesoe
        if(!ancestor_in_BAM){
            cerr << "Error: No data for ancestral sample '" << anc_tag <<
                    "' in the specifified BAM file. Check the sample tags match" << endl;
            exit(1);
        }
        if(!keepers.size()==0){
            cerr << "Sample(s) note persent in BAM file: ";
            for(auto s: keepers){
                cerr << s;
            }
            cerr << endl;
            exit(1)                
        }
        // And now.. go back over the read groups to map RG->sample index
        SampleMap samples;
        for(auto it = header.ReadGroups.Begin(); it!= header.ReadGroups.End(); it++){
            if(it->HasSample()){
                samples[it->ID] = name_map[it->Sample];
            }
        }
        return samples;

    }


    void ValidateNfreqs(boost::any& v, const vector<string>& values,  vector<double>*){
        vector<double> nfreqs {0,0,0,0};
        for(auto it: values){
            stringstream ss(*it);
            copy(istream_iterator<double>(ss), istream_iterator<double>(),
            back_inserter(nfreqs));        
        }
        if(nfeqs.size() != 4){
          throw po::validation_error("Must specify 4 (and only 4) nucleotide frequencies");        
        }
        v = nfreqs;
    }
    
    void check_args(po::variable_map &vm){
        // Is the experimental design one of the ones we can handle?
        if (vm["ploidy-ancestor"].as<int>() > 2 or vm["ploidy-ancestor"].as<int>() < 1){
            throw po::validation_error("accuMUlate can't only deal with haploid or diploid ancestral samples"); 
        }
        if (vm["ploidy-descendant"].as<int>() > 2 or vm["ploidy-desendant"].as<int>() < 1){
            throw po::validation_error("accuMUlate can't only deal with haploid or descendant samples");        
        }
        if (vm["ploidy-ancestor"].as<int>() == 1 and vm["ploidy-descendant"].as<int>() == 2){
            throw po::validation_error("accuMUlate has no model for a haploid->diploid MA experiemt");
        }
        //Do we have the right over-dispersion params set
        if (vm["ploidy-ancestor"].as<int>() == 1 or vm["ploidy-descendant"].as<int>() == 1){
            if(not vm.count("phi-haploid")){
                throw po::validation_error("Must specify phi-haploid (overdispersion for haploid sequencing)");
            }
        }
        if (vm["ploidy-ancestor"].as<int>() == 2 or vm["ploidy-descendant"].as<int>() == 2){
            if(not vm.count("phi-diploid")){
                throw po::validation_error("Must specify phi-diploid (overdispersion for diploid sequencing)");
            }
        }
    }


    void ParseCommandLineInput(int argc, char **argv, boost::program_options::variables_map &vm) {
//    boost::program_options::variables_map vm;
        namespace po = boost::program_options;
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
            ("out,o", po::value<string>()->default_value("acuMUlate_result.tsv"), "Out file name")
            ("intervals,i", po::value<string>(), "Path to bed file")
            ("config,c", po::value<string>(), "Path to config file")
            ("theta", po::value<double>()->required(), "theta")
            ("nfreqs", po::value<vector<double> >()->multitoken(), "Nucleotide frequencies")
            ("mu", po::value<double>()->required(), "Experiment-long mutation rate")
            ("seq-error", po::value<double>()->required(), "Probability of sequencing error")
            ("ploidy-ancestor", po::value<int>()->default_value(2), "Polidy of ancestor (1 or 2)")
            ("ploidy-descendant", po::value<int>()->default_value(2), "Ploidy of descendant (1 or 2)")
            ("phi-haploid",     po::value<double>(), "Over-dispersion for haploid sequencing")
            ("phi-diploid",     po::value<double>(), "Over-dispersion for diploid sequencing");
        po::variables_map vm;
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
            BamTools::SamHeader &header, BamTools::Fasta &reference_genome ) {


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
            vm["nfreqs"].as<vector< double> >(),
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
