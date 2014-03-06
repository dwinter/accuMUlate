sak - Slicing and dicing of FASTA/FASTQ files..
===============================================

SYNOPSIS
    sak [OPTIONS] [\B-o OUT.{fa,fq}] IN.{fa,fq}

DESCRIPTION
    "It slices, it dices and it makes the laundry!"

    Rewrite of the original SAK tool by Manuel Holtgrewe.

    -h, --help
          Displays this help message.
    --version
          Display version information

  Output Options:
    -o, --out-path FASTX
          Path to the resulting file. If omitted, result is printed to stdout.
          Use files ending in .fq or . to write out FASTQ. Valid filetypes
          are: .fq, .fastq, .fa, .fasta, .faa, .ffn, .fna, and .frn.
    -rc, --revcomp
          Reverse-complement output.
    -l, --max-length LEN
          Maximal number of sequence characters to write out.

  Filter Options:
    -s, --sequence NUM
          Select the given sequence for extraction by 0-based index.
    -sn, --sequence-name NAME
          Select sequence with name prefix being NAME.
    -ss, --sequences RANGE
          Select sequences from-to where from and to are 0-based indices.
    -i, --infix RANGE
          Select characters from-to where from and to are 0-based indices.
    -ll, --line-length LEN
          Set line length in output file. See section Line Length for details.
          In range [-1..inf].

LINE LENGTH
    You can use the setting --line-length for setting the resulting line
    length. By default, sequences in FASTA files are written with at most 70
    characters per line and sequences in FASTQ files are written without any
    line breaks. The quality sequence in FASTQ file is written in the same way
    as the residue sequence.

    The default is selected with a --line-length value of -1 and line breaks
    can be disabled with a value of 0.

USAGE EXAMPLES
    sak -s 10 IN.fa
          Cut out 11th sequence from IN.fa and write to stdout as FASTA.
    sak -ss 10-12 -ss 100-200 IN.fq
          Cut out 11th up to and including 12th and 101th up to and including
          199th sequence from IN.fq and write to stdout as FASTA.

VERSION
    sak version: 0.2
    Last update November 2012
