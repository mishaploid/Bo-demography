#!/usr/bin/env python3
# this bootstraps using the observed SFS 

import argparse
import json
import logging

from typing import Dict

import dadi

logging.basicConfig(level=logging.INFO)


def _cli() -> Dict[str, str]:
    """
    Define the command line arguments for building the SFS.

    Returns
    -------
    Key-value pairs with the arguments and their associated values passed
    at the command line as strings.
    """
    parser = argparse.ArgumentParser(
        description="Build the SFS from a VCF and pop-info file.",
        add_help=True,
    )
    required = parser.add_argument_group("required")
    required.add_argument(
        "-v",
        "--vcf",
        action="store",
        type=str,
        required=True,
        metavar="\b",
        help="Path to VCF file",
    )
    required.add_argument(
        "-p",
        "--pop_info",
        action="store",
        type=str,
        required=True,
        metavar="\b",
        help="Path to pop-info file",
    )
    required.add_argument(
        "-mc",
        "--model_config_file",
        action="store",
        type=str,
        required=True,
        metavar="\b",
        help="Path to model configuration json"
    )
    required.add_argument(
        "-m",
        "--model",
        action="store",
        type=str,
        required=True,
        metavar="\b",
        help="String for model specification"
    )
    additional = parser.add_argument_group("additional")
    additional.add_argument(
        "--suffix",
        action="store",
        type=str,
        metavar="\b",
        help="Additional suffix to add to file name",
    )
    additional.add_argument(
        "--subsample",
        action="store_true",
        help="Should the populations be subsampled?",
    )
    additional.add_argument(
        "--polarized",
        action="store_true",
        help="Should the SFS be polarized?",
    )
    args = vars(parser.parse_args())
    return args


def main() -> None:
    args = _cli()
    model_config = json.loads(open(args["model_config_file"]).read())
    sampling = (
        model_config["subsampling"] if args["subsample"] else model_config["sampling"]
    )
    logging.info("Reading in VCF file for bootstrapping...")
    bootstraps = dadi.Misc.bootstraps_subsample_vcf(args["vcf"], args["pop_info"], subsample = sampling, Nboot = 10, chunk_size = 5000000, pop_ids = model_config["models"][args["model"]]["subpops"], filter = True, flanking_info = [None, None], mask_corners = True, polarized = False)

    logging.info("Writing bootstrapped SFS files...")
    # i holds index for list length
    # rep holds actual fs object 
    for i,rep in enumerate(bootstraps): 
        rep.to_file(f"models/dadi/sfs_bootstrap/{args['model']}_{i}.fs")

if __name__ == "__main__":
    main()
