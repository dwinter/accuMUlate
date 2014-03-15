#Calling mutations from MA lines


##Required libs

Needs to be compiled/on path:
* [bamtools](https://github.com/pezmaster31/bamtools)
* [eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
* [Boost::program_options](http://www.boost.org/doc/libs/1_55_0/doc/html/program_options.html)

Provided in this repo
* bamtools headers  (for bamtools-utils, not in the normal include headers) 

##Compilation

At present `cmake` handles a basic build. Presuming bamtools and eigein are on
your path you can do this

```sh
cmake .
make
./run_tetma --config params.ini -b test/test.bam -r test/test.fasta -o test/test.out -p 0
```

The last line will  run the caller on a test dateset with ~2 000 bases, printing
the probability for each base

##TODO


* ~~have the main loop run the analysis from~~
    * ~~BAM file~~
    * ~~kist of sample IDs (one ancestor, rest descendants)~~
    * ~~Model Params:~~
        *~~theta~~
        *~~mu~~
        *~~nfreqs~~
        *~~phi-haploid~~
        *~~phi-diploid~~
    * ~~[BED file of coords to include]~~
    * ~~[probability threshold]~~
    * ~~[quality cut-off]~~
* ~~Boost PO/YAML/some other format to handle run variables~~
* ~~Manage build across platforms w/ CMake~~

* Optimize model code
* Clean up sample mapping/base calling
* Post-processor for putative mutations
* Create `.vcf` and annotations files from snps
* Make real test suite

