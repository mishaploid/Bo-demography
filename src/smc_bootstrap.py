#!/bin/python3

# script from: https://raw.githubusercontent.com/twooldridge/misc/master/SMC_bootstrap_BW.py

# This script is a modification of the SUPER HELPFUL script found in this github issue: https://github.com/popgenmethods/smcpp/issues/37.
# I've modified it to deal with irregularities in the input file, say due to some of the features of masking that vcf2smc produces


import click
import argparse
import os
import sys
import random
import gzip
import time


@click.command()
@click.option('--nr_bootstraps',  type=int, help="nr of bootstraps [20]", default=20)
@click.option("--chunk_size", type=int, help="size of bootstrap chunks [5000000]", default=5000000)
@click.option("--chunks_per_chromosome", type=int,
              help="nr of chunks to put on one chromosome in the bootstrap [20]", default=20)
@click.option("--nr_chromosomes", type=int, help="nr of chromosomes to write [30]", default=24)
@click.option("--seed", type=int, help="initialize the random number generator", default=None)
@click.argument("out_dir_prefix")
@click.argument("files", nargs=-1)
def main(nr_bootstraps, chunk_size, chunks_per_chromosome, nr_chromosomes, seed, out_dir_prefix, files):
    chunks = []
    offset = 0
    chunks_in_chrom = []
    if not seed:
        seed = int(time.time())
    random.seed(seed)
    print('seed: %s' % seed)
    for file in files:
        print(file)
        with gzip.open(file, 'rb') as f:
            header = f.readline().decode()
            prev_chunk_idx = -1
            for line_count,line in enumerate(f):
                line = line.decode()
                line = [int(x) for x in line.strip().split()]
                if line_count == 0:
                    pos = 1
                    offset += 1
                else:
                    pos = line[0] + offset
                    offset += line[0]
                chunk_index = (pos - 1) // chunk_size
                # The code below won't work if there are super large hom. stretches in the smc input file
                # This might happen if, for example you're only analyzing part of the chromosome, and you
                # mask everything outside of that part. The last line in the input file might look like:
                # 72402527 -1 0 0
                # Which tells SMC++ to mask the final ~72MB of the chromosome.
                # I've written a check below that should deal with this and stop the for loop if the next
                # position is more than (2*chunk_size) away
                if (chunk_index - prev_chunk_idx) >= 2:
                    print('\n')
                    print('Next position is in %s bases' % pos)
                    print('Last relative position was at %s bases' % lastpos)
                    print('This might be the result of how vcf2smc masks large stretches of the vcf (say for subsampling a chromosome), or it might be a real error.\nDouble check that this input file is formatted how you intend.\nStopping "chunking" of input file and moving on\n')
                    break
                print('Processing chunk %s of input file' % chunk_index)
                ##
                if chunk_index > len(chunks_in_chrom)-1:
                    chunks_in_chrom.append([])
                chunks_in_chrom[chunk_index].append(line)
                ##
                prev_chunk_idx = prev_chunk_idx + (chunk_index - prev_chunk_idx)
                lastpos = pos

        chunks.extend(chunks_in_chrom)

    for bootstrap_id in range(1, nr_bootstraps +1):
        for chr_ in range(1, nr_chromosomes + 1):
            chr_dir = "{}_{}".format(out_dir_prefix, bootstrap_id)
            if not os.path.exists(chr_dir):
                os.makedirs(chr_dir)
            chr_file = "{}/bootstrap_chr{}.gz".format(chr_dir, chr_)
            print("writing", chr_file, file=sys.stderr)
            with gzip.open(chr_file, 'wb') as f:
                f.write(header.encode())
                for i in range(chunks_per_chromosome):
                    chunk_id = random.randrange(len(chunks))
                    for line in chunks[chunk_id]:
                        line = ' '.join([str(x) for x in line]) + '\n'
                        f.write(line.encode())



if __name__ == "__main__":
    main()
