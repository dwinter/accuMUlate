#!/bin/sh
echo '===running accumulate on test data==='
echo '(bamtools is going to complain about the not-very-well sorted BAM)'
build/accuMUlate -c example/test_params.ini -b example/test.bam -r example/test.fasta
-o example/test.out 

echo '\n===First pass result ==='
cat test/test.out

echo '\n===running post processor==='
build/pp -i test/test.out -b test/test.bam -c test/pp_params.ini -o  test/filter.out

echo '\n===filter results==='
cat test/filter.out
