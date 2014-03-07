#Calling mutations from MA lines


##Required libs

Needs to be compiled/on path:
* [bamtools](https://github.com/pezmaster31/bamtools)
* [eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)

Provided in this repo
* bamtools headers        
* [SeqAn](https://www.seqan.de/)

##Compilation

At the moment somehting like this (presuming `eigen` and the compiled `bamtools`
libs are on your path):

`clang++ model.cc parse_bam.cc -lbamtools -std=c++11 -lbamtools-utils -Ithird-party/ -Ithird-party/bamtools/src -g -Wall -o run_tetma`

This throws around 40 warnings on my system (mostly related to `SeqAn`).

In the future: use CMake to handle build.

##TODO
* have the main loop run the analysis from
    * BAM file
    * List of sample IDs (one ancestor, rest descendants)
    * Model Params:
        theta
        mu
        nfreqs
        phi-haploid
        phi-diploid
    * [BED file of coords to include]
    * [probability threshold]
    * [quality cut-off]
* Boost PO/YAML/some other format to handle run variables
* Manage build across platforms w/ CMake
