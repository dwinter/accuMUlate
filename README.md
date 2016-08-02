[![Build Status](https://travis-ci.org/dwinter/accuMUlate.svg?branch=master)](https://travis-ci.org/dwinter/accuMUlate)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.19942.svg)](http://dx.doi.org/10.5281/zenodo.19942)
#Calling mutations from MA lines

`accuMUlate` is a mutation caller designed with Mutation Accumulation (MA)
experiments in mind. Up until recently the software has been developed with the
goal of studying [a specific MA experiment](http://biorxiv.org/content/early/2015/09/10/025536),
 but we are in the process of adapting the program to work on any MA dataset you might have. 

Please [contact us](mailto:david.winter@gmail.com) if you are interested in using the software
while we develop it.


##Required libs

Needs to be compiled/on path:
* [bamtools](https://github.com/pezmaster31/bamtools)
* [eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
* [Boost::program_options](http://www.boost.org/doc/libs/1_55_0/doc/html/program_options.html)

~~Provided in this repo~~
~~* bamtools headers  (for bamtools-utils, not in the normal include headers)~~

##Compilation

At present `cmake` handles a basic build. Presuming bamtools and eigein are on
your path you can do this

```sh
cmake .
make
./accuMUlate -c example/test_params.ini -b example/test.bam -r example/test.fasta -o example/test.out
./pp -c example/pp_params.ini -b example/test.bam -i example/test.out 

```

The last line will  run the caller on a test dateset with 6 000 bases, and
show find mutations in each gene.

