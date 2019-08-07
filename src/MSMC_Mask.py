import gzip
import argparse


def main():
    # Parse input arguments
    parser = argparse.ArgumentParser(description="script to generate the mask file from the vcf file generated with GATK for MSMC")
    parser.add_argument("-v", "--vcf", action="store", required=True, help="Input VCF file. Should be a unisample-unichromosome vcf")
    parser.add_argument("-o", "--out", action="store", required=True, help="Output filename")
    parser.add_argument("-d", "--depth", type=float, action="store", required=True, help="Median depth")
    parser.add_argument("-g", "--gzip", action="store_true", required=False, help="Set if the VCF is gzipped.")

    args = parser.parse_args()

    vcf_in = args.vcf
    out_name = args.out
    depth = args.depth

    # Open the input vcf file
    if args.gzip:
        opener = gzip.open
    else:
        opener = open

    f_out = open(out_name, 'w')
    minDepth = depth / 2.0
    maxDepth = depth * 2.0

    # Read the input vcf file
    with opener(vcf_in, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                pass
            else:
                # Parse relevant vcf columns
                line = line.split("\t")
                chrom = line[0]
                pos = int(line[1])
                info = line[9:]
                infos = info[0].split(":")
                if infos[2] != ".":
                    dp = int(infos[2])
                else:
                    dp = 0
                if dp >= minDepth and dp <= maxDepth:
                    p2=pos+1
                    f_out.write("{}\t{}\t{}\n".format(chrom, pos, p2))
                else:
                    pass
if __name__ == "__main__":
    main()
