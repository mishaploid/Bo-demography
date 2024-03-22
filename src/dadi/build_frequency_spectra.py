#!/usr/bin/env python3

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
    model_config = json.loads(open("src/dadi/model_config.json").read())
    sampling = (
        model_config["subsampling"] if args["subsample"] else model_config["sampling"]
    )
    logging.info("Reading in VCF file...")
    if args["subsample"]:
        data_dict = dadi.Misc.make_data_dict_vcf(
            args["vcf"], args["pop_info"], subsample=sampling
        )
    else:
        data_dict = dadi.Misc.make_data_dict_vcf(args["vcf"], args["pop_info"])

    logging.info("Building frequency spectra for each model...")
    for fs in model_config["models"]:
        model_sub = model_config["models"][fs]
        print(f"Building sfs for populations: {', '.join(model_sub['subpops'])}")
        sfs = dadi.Spectrum.from_data_dict(
            data_dict=data_dict,
            pop_ids=model_sub["subpops"],
            projections=[2 * sampling[pop] for pop in model_sub["subpops"]],
            polarized=args["polarized"]
        )
        sfs.to_file(f"models/dadi/sfs/{model_sub['sfs_file']}")


if __name__ == "__main__":
    main()
