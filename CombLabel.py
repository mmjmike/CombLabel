#!/usr/bin/python3

from classes.ucsl_io import find_ncs, read_prices, read_sequence, read_stock, read_assignment
from classes.logger import create_logger_main
from classes.search_objects import CLSequence, Stock
from classes.sequence_specific import CLOptimizer
import argparse
import os


DEFAULT_LOG_FILENAME = "CombLabel.log"

script_path = os.path.split(os.path.realpath(__file__))[0]


def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("ncs", help='Specify NCS name or filename', type=str)
    parser.add_argument("sequence_file", help='Specify sequence filename', type=str)
    parser.add_argument("stock_file", help='Specify amino acid stock filename', type=str)
    parser.add_argument("--price", "-p", help='Specify price file', default="", type=str)
    parser.add_argument("--jobname", "-j",
                        help="Jobname for this calculation", default="",
                        type=str)
    parser.add_argument("--number_samples", "-n",
                        help="Specify starting number of samples (1-9)", default=1,
                        type=int)
    parser.add_argument("--assignment", "-a",
                        help="Specify the filename containing the numbers of already assigned residues", default="",
                        type=str)
    parser.add_argument('--verbose', '-v',
                        help='Verbose output to console',
                        action="store_true")
    parser.add_argument('--silent', '-s',
                        help='No output to console (overrides verbosity)',
                        action="store_true")
    return parser.parse_args()


def read_parameters(args, logger):
    logger.info("Reading NCS")
    ncs, msg = find_ncs(args.ncs, script_path)
    if not ncs:
        logger.error(msg)
        exit()
    logger.info("NCS: {}".format(ncs.name))

    sequence = ""
    logger.info("Reading sequence: '{}'".format(args.sequence_file))
    if os.path.isfile(args.sequence_file):
        sequence, msg = read_sequence(args.sequence_file)
        if not sequence:
            logger.error(msg)
            exit()
    else:
        logger.error("Error! Sequence file '{}' not found".format(args.sequence_file))

    if len(sequence) < 3:
        print("Sequence length ({}) is too short".format(len(sequence)))
        exit()
    logger.info("Sequence length - {} residues".format(len(sequence)))


    stock = {}
    if os.path.isfile(args.stock_file):
        logger.info("Reading stock: '{}'".format(args.stock_file))
        stock, msg = read_stock(args.stock_file)
        if not stock:
            logger.error(msg)
            exit()
    else:
        logger.error("Error! Stock file '{}' not found".format(args.stock_file))
    logger.info("Stock - OK")


    optimize_price = True
    if not args.price:
        optimize_price = False
    prices = {}
    if optimize_price:
        logger.info("Reading prices: '{}'".format(args.price))
        if os.path.isfile(args.price):
            prices, msg = read_prices(args.price)
            if not prices:
                logger.error(msg)
                exit()
        else:
            logger.error("Error! Price file '{}' not found".format(args.price))
            exit()
        logger.info("Prices - OK")


    assignment = True
    if not args.assignment:
        assignment = False
    assignment_numbers = set()
    if assignment:
        logger.info("Reading assignment: '{}'".format(args.assignment))
        if os.path.isfile(args.assignment):
            assignment_numbers, msg = read_assignment(args.assignment)
        else:
            logger.error("Error! Assignment file '{}' not found".format(args.assignment))
    for assignment_num in assignment_numbers:
        if assignment_num > len(sequence):
            warning_msg = "Warning! Residue number {} is exceeding sequence length ({})".format(assignment_num,
                                                                                                len(sequence))
            logger.error(warning_msg)
    if assignment:
        logger.info("{} backbone assignments read".format(len(assignment_numbers)))

    jobname = args.jobname
    seq_name = args.sequence_file.split(".")[0]
    if not jobname:
        jobname = "{}_{}".format(seq_name, ncs.name)


    samples = args.number_samples # default value is 1
    if samples < 1 or samples > 9:
        msg = "Error! Number of samples is not in range 1-9"
        logger.error(msg)
        exit()

    parameters = {
        "ncs": ncs,
        "sequence": sequence,
        "stock": stock,
        "optimize_price": optimize_price,
        "prices": prices,
        "samples": samples,
        "jobname": jobname,
        "verbose": args.verbose,
        "silent": args.silent,
        "assignment": assignment_numbers
    }
    return parameters


def find_solution(parameters, logger):
    sequence = CLSequence(parameters["sequence"])
    stock = Stock(parameters["stock"], parameters["prices"])
    sequence.sequence_stats(stock, assignment=parameters["assignment"])
    scheme_optimizer = CLOptimizer(parameters, logger, sequence)
    scheme_optimizer.run()
    # scheme_optimizer.write_results()


def main():
    args = read_args()
    logger = create_logger_main(args, DEFAULT_LOG_FILENAME)
    parameters = read_parameters(args, logger)
    find_solution(parameters, logger)


if __name__ == "__main__":
    main()
