#!/usr/bin/python3

import argparse
import os
import logging
from classes.search_objects import BlockFinder
from classes.ucsl_io import find_ncs, write_ncs_stamp
from classes.interactive import answer_yes
import re
import cProfile

# from classes.interactive import ask_to_continue_input

script_path = os.path.split(os.path.realpath(__file__))[0]


DEFAULT_LOGFILE_NAME = "blockfinder.log"


def good_filename(filename):
    return len(filename)


def run_blockfinder_once(parameters, logger, elb_logger):
    ncs = parameters["ncs"]
    samples = parameters["samples"]
    min_depth = parameters["min_patterns"]
    min_t_free_value = parameters["min_t_free"]

    elb_logger.info(write_ncs_stamp(ncs))
    all_blocks = {}
    block_finder = BlockFinder(samples, ncs, min_depth,  logger, elb_logger,
                               block_finder_mode=True, min_t_free = min_t_free_value)
    block_finder.find()
    result = block_finder.result
    all_blocks.update({samples: result})


def find_blocks(parameters, logger, elb_logger):
    ncs = parameters["ncs"]
    max_block_size = parameters["samples"]
    elb_logger.info(write_ncs_stamp(ncs))
    # self.output("Search for elementary blocks with max {} samples".format(self.max_block_size))
    min_depth = 1
    all_blocks = {}
    for i in range(max_block_size):
        # self.output("Search in {} samples:".format(i+1))
        # self.output("Time in find_blocks: {} seconds\n".format(str(time.time()-t0)))
        block_finder = BlockFinder(i + 1, ncs, min_depth,  logger, elb_logger)
        block_finder.find()
        result = block_finder.result
        all_blocks.update({i+1: result})
        if result:
            min_depth = max(result) + 1
        else:
            min_depth += 1


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("ncs", help='Specify NCS (NMR Coding System)', type=str)
    parser.add_argument("samples",
                        help='Specify maximal number of samples (range 1-9) or exact number of samples (if run in scheme finding mode)',
                        type=int)
    parser.add_argument("--minpatterns", "-p",
                        help='Specify minimal number of patterns to run in scheme finding mode',
                        type=int)
    parser.add_argument('--verbose', '-v',
                        help='Verbose output to console',
                        action="store_true")
    parser.add_argument('--silent', '-s',
                        help='No output to console (overrides verbosity)',
                        action="store_true")
    parser.add_argument("--mintfree", "-t", dest="min_t_free",
                        help='Required number of patterns without \"T\" labeling type in ELB', default=-1, type=int)
    return parser.parse_args()


def get_params(args, logger):
    ncs, msg = find_ncs(args.ncs, script_path)
    if not ncs:
        logger.error(msg)
        exit()
    else:
        logger.info(msg)

    samples = args.samples
    if samples < 1 or samples > 9:
        logger.error("Number of samples '{}' is not in range 1-9".format(samples))
        exit()

    min_patterns = 0
    exact_patterns = False

    if args.minpatterns:
        min_patterns = args.minpatterns
        exact_patterns = True

    min_t_free = -1
    check_t_free = False
    if args.min_t_free:
        if args.min_t_free >=0:
            min_t_free = args.min_t_free
            check_t_free = True
        else:
            min_t_free = -1

    if check_t_free:
        exact_patterns = True

    params = {
        "ncs": ncs,
        "samples": samples,
        "exact_patterns": exact_patterns,
        "min_patterns": min_patterns,
        "check_t_free": check_t_free,
        "min_t_free": min_t_free,
        "verbose": args.verbose,
        "silent": args.silent
    }
    return params


def make_loggers(args):
    logger = logging.getLogger("Logfile_logger")
    logger.setLevel(logging.DEBUG)

    min_patterns = 0
    exact_patterns = False

    if args.minpatterns:
        min_patterns = args.minpatterns
        exact_patterns = True

    filename = args.ncs
    filename += "_" + str(args.samples)
    if exact_patterns:
        filename += "_" + str(min_patterns)
    log_filename = "{}_blocks.log".format(filename)
    if os.path.isfile(log_filename):
        input_key = input("Logfile '{}' already exists."
                          "Do you want to overwrite? (Y[es]/N[o]) > ".format(log_filename))
        if not answer_yes(input_key):
            exit()

    logfile_handler = logging.FileHandler(log_filename, mode='w')
    logfile_format = logging.Formatter('%(asctime)s:%(message)s')
    logfile_handler.setFormatter(logfile_format)
    logfile_handler.setLevel(logging.DEBUG)
    stream_handler = logging.StreamHandler()
    stream_level = logging.INFO
    if args.verbose:
        stream_level = logging.DEBUG
    if args.silent:
        stream_level = logging.CRITICAL
    stream_handler.setLevel(stream_level)
    logger.addHandler(logfile_handler)
    logger.addHandler(stream_handler)


    elb_filename = "{}.elb".format(filename)
    elb_logger = logging.getLogger("ELB_Logger")
    elb_file_handler = logging.FileHandler(elb_filename, mode='w')
    elb_logger.setLevel(logging.INFO)
    elb_logger.addHandler(elb_file_handler)


    return logger, elb_logger


def main():

    args = get_args()
    logger, elb_logger = make_loggers(args)
    parameters = get_params(args, logger)
    if parameters["exact_patterns"]:
        run_blockfinder_once(parameters, logger, elb_logger)
    else:
        find_blocks(parameters, logger, elb_logger)

    print("Finished successfully")


if __name__ == "__main__":
    main()


    # while True:
    #     input_filename = input("Please input filename for BlockFinder results: ")
    #     if not good_filename(input_filename):
    #         print("'{}' can't be used as a filename. Please enter a valid one (0 < length < 256; can't use following symbols: )")
    #         if not ask_to_continue_input():
    #             return
