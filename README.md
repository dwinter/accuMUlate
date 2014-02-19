Software for calling mutations from MA experiments in Thermophillia

#The plan:

* move Reed's model/calc code into a seperate file to include
* have the main loop run the analysis from
    ** BAM file
    ** List of sample IDs (one ancestor, rest descendants)
    ** Model Params:
        theta
        mu
        nfreqs
        phi-haploid
        phi-diploid
    ** [BED file of coords to include]
    ** [probability threshold]
    ** [quality cut-off]
* Bost/YAML/some sort of parser for these inputs
* Manage build w/ CMake
