#Calling mutations from MA lines


##Required libs

Needs to be compiled/on path:
* [bamtools](https://github.com/pezmaster31/bamtools)
* [eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)

Provided in this repo
* bamtools headers        
* [SeqAn](https://www.seqan.de/)

##Compilation

At present `cmake` handles a basic build. Presuming bamtools and eigein are on
your path you can do this

```sh
cmake .
make
./run_tetma
```
The buld process will throw a lot of warnigns (mainly form `seqan`). The last
run the caller on a test dateste with ~2 000 bases.

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
