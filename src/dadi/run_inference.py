#!/usr/bin/env python3

import argparse
import json
import logging
from typing import Dict

import dadi
import dadi.NLopt_mod

from demographic_models import models

logging.basicConfig(level=logging.INFO)


def _cli() -> Dict[str, str]:
    """
    Define the command line arguments for running demographic inference.

    Returns
    -------
    Key-value pairs with the arguments and their associated values passed
    at the command line as strings.
    """
    parser = argparse.ArgumentParser(
        description="Run inference for a demographic model",
        add_help=True,
    )
    required = parser.add_argument_group("required")
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
        "--model",
        action="store",
        type=str,
        required=True,
        metavar="\b",
        help="Name of model to run",
    )

    args = vars(parser.parse_args())
    return args


def main() -> None:
    args = _cli()
    model_config = json.loads(open(args["model_config_file"]).read())
    model_names = [m for m in model_config["models"]]
    
    # Check to make sure that 'model' argument
    # is in the list of models
    if args["model"] not in model_names:
        print(f"Could not find model named {args['model']}...")
        print(
            f"Please choose a model from the following list:\n{', '.join(model_names)}"
        )
        raise RuntimeError

    model_info = model_config["models"][args["model"]]

    data = dadi.Spectrum.from_file(f"models/dadi/sfs/{model_info['sfs_file']}")
    ns = data.sample_sizes
    pts_l = [max(ns) + i for i in [10, 20, 30]]

    upper_bound = model_info["upper_bound"]
    lower_bound = model_info["lower_bound"]

    # Initial guess for parameter values (will be randomly perturbed)
    p_init = model_info["p_init"]
    model_func = models[args["model"]]

    with open(f"models/dadi/results/{args['model']}.csv", "w") as f_out:
        for r in range(model_config["optimization_reps"]):
            if r % 10 == 0:
                logging.info(f"\n\n*** Starting rep {r + 1} ***\n\n")

            p0 = dadi.Misc.perturb_params(
                p_init, fold=2, lower_bound=lower_bound, upper_bound=upper_bound
            )

            popt, LLopt = dadi.Inference.opt(
                p0,
                data,
                model_func,
                pts_l,
                lower_bound=lower_bound,
                upper_bound=upper_bound,
                verbose=20,
            )

            model = model_func(popt, ns, pts_l)
            theta = dadi.Inference.optimal_sfs_scaling(model, data)
            line_buffer = [str(x) for x in ([r + 1, LLopt] + list(popt) + [theta])]
            print(",".join(line_buffer), file=f_out)


if __name__ == "__main__":
    main()
