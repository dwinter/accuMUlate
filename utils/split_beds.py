#! /usr/bin/env python

"""
Usage 

split_beds [bam] [# of files to split into] [directory to write bed files]

requires samtools to be on $PATH... everyone has samtools right?

"""

import sys
import subprocess
import collections

class Chromosome(object):
    """ Represents one chromosome from a BAM header """

    def __init__(self, name, length, end):
        self.name = name
        self.length = length
        self.end = end

    def make_bed_line(self, bed_start=1, bed_end=None):
        if not bed_end:
            bed_end = self.length
        return("{0}\t{1}\t{2}\n".format(self.name, bed_start, bed_end))


def parse_header(fname):
    cmd = "samtools view -H {0}".format(fname)
    child = subprocess.Popen(cmd, stdout=subprocess.PIPE, 
                                  stderr=subprocess.PIPE,
                                  shell=(sys.platform!="win32"))
    header, err = child.communicate()
    genome = []
    cummulative_len = 0
    for line in header.split("\n"):
        if line.startswith("@SQ"):
            chrom, l = [e.split(":")[1] for e in line.split()[1:] ]
            cummulative_len += int(l)
            genome.append(Chromosome(chrom, int(l), cummulative_len))
    return(genome)
    

def main():
    genome = parse_header(sys.argv[1])
    out_dir = sys.argv[3]
    genome_len = genome[-1].end
    step = genome_len/int(sys.argv[2])
    genome_len = genome[-1].end
    bed_ends = collections.deque(range(step, genome_len, step) + [genome_len, None])

    start = 1
    chr_idx = 0
    e = bed_ends.popleft()
    fcounter = 0

    while bed_ends:
        fcounter += 1
        with open("{0}/interval_{1}.bed".format(out_dir, fcounter), "w") as out:
            while True:
                if e > genome[chr_idx].end:
                    out.write(genome[chr_idx].make_bed_line(bed_start=start))
                    chr_idx += 1
                    start =  1
                else:
                    offset = e - (genome[chr_idx].end - genome[chr_idx].length)
                    out.write(genome[chr_idx].make_bed_line(bed_start = start,
                        bed_end = offset))
                    start = offset + 1
                    e = bed_ends.popleft()
                    break

    print("wrote {0} bed files for a total of {1} bases".format(fcounter, genome_len) )
    return 0

if __name__ == "__main__":
    main()


