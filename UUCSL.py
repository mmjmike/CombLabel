from classes.search_objects import SchemeOptimizer
from classes.ucsl_io import read_blocks, find_ncs, read_prices
import argparse
import os
import logging

script_path = os.path.split(os.path.realpath(__file__))[0]


def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("ncs", help='Specify NCS name or filename', type=str)
    parser.add_argument("elb_file", help='Specify ELB file', type=str)
    parser.add_argument("price_file", help='Specify price file', type=str)
    parser.add_argument("--jobname", "-j",
                        help="Jobname for this calculation", default="",
                        type=str)
    parser.add_argument('--verbose', '-v',
                        help='Verbose output to console',
                        action="store_true")
    parser.add_argument('--silent', '-s',
                        help='No output to console (overrides verbosity)',
                        action="store_true")
    return parser.parse_args()


def check_args():
    pass


def read_parameters(args):
    ncs = find_ncs(args.ncs, script_path)
    result, blocks, ncs_name = read_blocks(args.elb_file)
    if result:
        print(result)
        exit()
    prices = read_prices(args.price_file)
    max_block_samples = max(blocks.keys())
    jobname = args.jobname
    if not jobname:
        jobname = "{}_{}".format(ncs.name, max_block_samples)
    parameters = {
        "blocks": blocks,
        "ncs": ncs,
        "prices": prices,
        "jobname": jobname,
        "verbose": args.verbose,
        "silent": args.silent
    }
    return parameters


def optimize_scheme(parameters, logger):
    scheme_optimizer = SchemeOptimizer(parameters, logger)
    scheme_optimizer.run()


def create_logger(parameters):
    logger = logging.getLogger("logger")
    return logger

def main():
    args = read_args()
    parameters = read_parameters(args)
    logger = create_logger(parameters)
    best_scheme = optimize_scheme(parameters, logger)
    write_best_scheme(best_scheme)


if __name__ == "__main__":
    main()