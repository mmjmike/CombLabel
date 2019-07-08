#!/usr/bin/python3

from classes.search_objects import SchemeOptimizer
from classes.ucsl_io import read_blocks, find_ncs, read_prices
from classes.constants import Constants
from classes.logger import create_logger_main
import argparse
import os


DEFAULT_LOG_FILENAME = "UUCSL.log"

script_path = os.path.split(os.path.realpath(__file__))[0]


def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("ncs", help='Specify NCS name or filename', type=str)
    parser.add_argument("elb_file", help='Specify ELB file', type=str)
    parser.add_argument("price_file", help='Specify price file', type=str)
    parser.add_argument("--aminoacids", '-a', help='Specify the set of amino acids to be used in the labeling scheme (one letter code without seperation, default is all 20 AA)', type=str)
    parser.add_argument("--notaminoacids", '-n',
                        help='Specify the set of amino acids that will not be used in the labeling scheme (one letter code without seperation)',
                        type=str)
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


def extract_residues(line):
    residues = []
    for char in line:
        res = char.upper()
        if res in residues:
            print("Warning! Residue '{}' is repeated at least twice in arguments".format(res))
        elif res in Constants.RES_TYPES_LIST:
            residues.append(res)
        else:
            print("Warning! Residue '{}' is unknown (And will be skipped in calculation0".format(res))
    return residues


def read_parameters(args, logger):
    ncs, msg = find_ncs(args.ncs, script_path)
    if not ncs:
        logger.error(msg)
        exit()
    result, blocks, ncs_name, deuterated = read_blocks(args.elb_file, logger)
    if result:
        exit()
    prices = {}
    if os.path.isfile(args.price_file):
        prices, msg = read_prices(args.price_file)
        if not prices:
            logger.error(msg)
            exit()
    else:
        logger.error("Error! Price file '{}' not found".format(args.price_file))
    max_block_samples = max(blocks.keys())
    jobname = args.jobname
    if args.aminoacids:
        residues = extract_residues(args.aminoacids)
    else:
        residues = Constants.RES_TYPES_LIST

    if args.notaminoacids:
        not_residues = extract_residues(args.notaminoacids)
    else:
        not_residues = []
    final_residues = []
    for res in residues:
        if res not in not_residues:
            final_residues.append(res)
    if not jobname:
        jobname = "{}_{}".format(ncs.name, max_block_samples)
    parameters = {
        "blocks": blocks,
        "ncs": ncs,
        "prices": prices,
        "residues": final_residues,
        "jobname": jobname,
        "verbose": args.verbose,
        "silent": args.silent
    }
    return parameters


def optimize_scheme(parameters, logger):
    scheme_optimizer = SchemeOptimizer(parameters, logger)
    scheme_optimizer.run()


def main():
    args = read_args()
    logger = create_logger_main(args, DEFAULT_LOG_FILENAME)
    parameters = read_parameters(args, logger)
    optimize_scheme(parameters, logger)


if __name__ == "__main__":
    main()
