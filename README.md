[![Build Status](https://travis-ci.org/dwinter/accuMUlate.svg?branch=master)](https://travis-ci.org/dwinter/accuMUlate)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.19942.svg)](http://dx.doi.org/10.5281/zenodo.19942)
#Calling mutations from MA lines

`accuMUlate` is a mutation caller designed with Mutation Accumulation (MA)
experiments in mind. The probablistic approach to mutation calling implemented by 
`accuMUlate` allows us to identify putative mutations, accomidate noise produced
by NGS sequencing and accomidate diploid, haploid or diploid to haploid
experimental designs.


##Getting started

The [wiki](https://github.com/dwinter/accuMUlate/wiki) gives a detailed account
of how to prepare your rate, compile accuMUlate, run the program and understand
results. If you want to get started even more quickly here's what you need to
know

##Prerequisites 

In order to install `accuMUlate` you will need the following libraries

* [bamtools](https://github.com/pezmaster31/bamtools)
* [eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
* [Boost::program_options](http://www.boost.org/doc/libs/1_55_0/doc/html/program_options.html)

You will also need [`CMake`](https://cmake.org/) to manage the build. Package
managers for linux distributions and OSX  should let you installed pre-compiled
version of Eigen, Boost and CMake. The [wiki](https://github.com/dwinter/accuMUlate/wiki) 
describes how to compile bamtools.

###Compilation

With prerequsites installed, building the software is easy:

```sh
cd build
cmake ..
make
```

###Test run

That above commands should make two programs in the `build` directory:
`accuMUlate` (the mutation called) and `denominate` (a tool for calculating the
number of callable sites in a BAM file). To test the everything has gone well
you can run these on some test data. First `accuMUlate`, which should produce
a warning message and detailed information about one possible mutation:


```sh
./accuMUlate -c test/data/example_params.ini \
             -b test/data/test.bam \
             -r test/data/test.fasta \
             -i test/data/test.bed 
```

```
Warning: excluding data from 'D6' which is included in the BAM file but not the list of included samples
good_mutation	600	601	C	D1	CC->G	0.999999	0.999999	1	1.86137e-10	137	13	11	0	011.8531	-0.756908	1	1	
```

And then denominate, the mysterious string of integer are teh number of
ancestrally "A", "C", "G" and "T" bases for each sample that could be called for
a mutation is one was present.

```sh
./denominate -c test/data/example_params.ini \
            -b test/data/test.bam \
            -r test/data/test.fasta \
            -i test/data/test.bed \
            --max-depth 150
```

```
Warning: excluding data from 'D6' which is included in the BAM file but not the list of included samples
2	1	0	0	2	1	0	0	2	1	0	0	2	1	0	0	2	1	0   0	
```

## How does it work

The `accuMUalte` model is descirbed in. 

Winter, Long et al (in review) Low base-substitution mutation rate in the ciliate _Tetrahymena thermophila_
pre-print doi:[https://dx.doi.org/10.1101/025536](https://dx.doi.org/10.1101/025536)

   

## Help, Bugs  and Suggestions.

Your first stop should be the [wiki](https://github.com/dwinter/accuMUlate/wiki), 
which contains information about the input and output files and more detail
about using `accMUlate`. If that doesn't help please file issues at this
repository or [email David Winter](mailto:david.winter@gmail.com).
