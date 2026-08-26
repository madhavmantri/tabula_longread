#!/usr/bin/env python

"""
Subsample sorted, corrected bam file to flatten coverage across cell barcodes.

Requires the python library pysam and samtools to be in your path.

"""

import argparse
import pysam
from collections import defaultdict
import random
import subprocess
import time
import os
import sys


# Command line options
parser = argparse.ArgumentParser(description='Subsample a sorted, corrected BAM file '
                                 'to cap the number of reads per cell barcode.')
parser.add_argument("-b","--bam",
					type=str,
					help="Input bamfile with cell barcode calling")
parser.add_argument("-m","--maxreads",
					type=int,
                    default=300000,
					help="maximum number of reads to retain per cell barcode")
parser.add_argument("-o","--output",
					type=str,
                    help="Path to output file.")
parser.add_argument("-t","--tag",
					type=str,
                    default="CB",
                    help="BAM tag for cell barcode (default: CB)")
parser.add_argument("-@","--threads",
					type=int,
                    default=1,
                    help="Number of threads for samtools (default: 1)")
args=parser.parse_args()

#print(args.description)

bamfile = args.bam
max_reads = args.maxreads
output = args.output
tag = args.tag
threads = args.threads

# ===                Functions                 ===

def time_log(func):
    """Decorator that measures how long a function takes to run."""
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        end = time.perf_counter()
        print(f"{func.__name__} took {end - start:.4f} seconds.")
        return result
    return wrapper

@time_log
def pull_info(bam, tag):
    """
    Iterates over reads to pull the read names and group them into
    a cell barcodes dictionary.

    Returns the dictionary with cell barcodes as keys and a list
    of reads with that cell barcode tag as values.
    """
    reads_tag = defaultdict(list)
    for read in bam.fetch(until_eof=True):
        if read.has_tag(tag):
            reads_tag[read.get_tag(tag)].append(read.query_name)
    return(reads_tag)

@time_log
def subsample_cbreads(reads_tag, max_reads):
    """
    Iterates over cell barcodes dictionary.

    If the number of reads
    associated with a given cell barcode is greater than the
    max_reads value, then a random subset of N=max_reads is selected.

    Returns flattened list of reads to retain for the subsample.
    """
    reads = []
    for t in reads_tag:
        if len(reads_tag[t]) > max_reads:
            sampled = random.sample(reads_tag[t], max_reads)
            reads.append(sampled)
        else:
            reads.append(reads_tag[t])
    reads = [cb_read for cb_reads in reads for cb_read in cb_reads]
    return(reads)


# ===                Execution                 ===
print("Starting execution")
bam = pysam.AlignmentFile(bamfile, "rb", check_sq=False)

print("Reading bam tags")
reads_tag = pull_info(bam, tag)

print("Selecting a random subsample for cell barcodes with over ",
      str(max_reads), " reads")
keep_reads = subsample_cbreads(reads_tag, max_reads)

tmpfile = os.path.join(os.path.dirname(output), "subsampled_reads.txt")
print("Writing subsampled read names to", tmpfile)
with open(tmpfile, "w") as file:
    for read in keep_reads:
        file.write(read + "\n")

print("Subsampling bam file with samtools")
command = f"samtools view -b -@ {threads} -o {output} --qname-file {tmpfile} {bamfile}"
ret = subprocess.call(command, shell=True)

os.remove(tmpfile)

if ret != 0:
    print("samtools view failed", file=sys.stderr)
    sys.exit(ret)
