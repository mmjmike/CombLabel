#!/usr/bin/python3 -u

import argparse
import os

from UUCSL.elbclean import clear_empty_block_types
from classes.constants import Constants, Pattern
from classes.ucsl_io import make_block_stats, add_to_file_name, write_blocks, read_blocks

LOG_ITERATION = 100000
index_of_typeT = Pattern.index_of_type(Constants.typeT)

def filterT_blocks(blocks, min_t_free, verbose_output):
    iter = 0
    total_blocks = make_block_stats(blocks)["total"]
    filtered_blocks = 0
    print("Filter blocks with less then {} patterns that have no \"T\" labeling".format(min_t_free))
    for samples_num in blocks:
        patterns_numbers = list(blocks[samples_num].keys())
        patterns_numbers.sort(reverse=True)
        for pattern_num in patterns_numbers:
            new_block_list = []
            for block in blocks[samples_num][pattern_num]:
                iter += 1
                if iter % LOG_ITERATION == 0:
                    print("filterT blocks checked {}/{}".format(iter, total_blocks))
                have_t, have_no_t = count_typeT(block)
                output = "have_t = {:2} have_no_t = {:2}, ".format(have_t, have_no_t)
                if have_no_t >= min_t_free:
                    output += "APPEND block"
                    new_block_list.append(block)
                    filtered_blocks += 1
                else:
                    output += "SKIP block"
                if verbose_output:
                    print(block.full_str())
                    print(output + "\n")
            if not new_block_list:
                blocks[samples_num].pop(pattern_num)
            else:
                blocks[samples_num][pattern_num] = new_block_list
    print("Filtered {}/{} blocks by T labeling, min_t_free = {}".format(filtered_blocks, total_blocks, min_t_free))
    return clear_empty_block_types(blocks)

def count_typeT(block):
    return Pattern.count_type_in_list_of_simplified(block.simplified, index_of_typeT)

def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("elb_files", help='Specify ELB file(s)', type=str, nargs='+')
    parser.add_argument("-o", dest="output_file", help='Specify output file (optional)', default="", type=str)
    parser.add_argument("--verbose", "-v", dest="verbose", help='verbose output', action = "store_true")
    parser.add_argument("--mintfree", "-t", dest = "min_t_free",
                        help='Required number of patterns without \"T\" labeling type in ELB', default=16, type=int)
    args = parser.parse_args()
    for file in args.elb_files:
        if not os.path.isfile(file):
            print("Error! ELB file '{}' not found".format(file))
            exit()
        print("ELB file: {}".format(file))
    if args.output_file:
        print("Output file: {}".format(args.output_file))
        output_file = args.output_file
    else:
        output_file = add_to_file_name(args.elb_files[0] + "_" + args.elb_files[-1], "_filtered")
    return args.elb_files, output_file, args.min_t_free, args.verbose

def main():
    elb_files, output_file, min_t_free, verbose_output = read_args()
    blocks = {}
    for file in elb_files:
        result, blocks, ncs_name, deuterated = read_blocks(file, initial_blocks = blocks)
        print("{} read".format(file))
    new_blocks = filterT_blocks(blocks, min_t_free, verbose_output)
    write_blocks(new_blocks, ncs_name, output_file, deuterated)
    print("Blocks written to '{}' successfully".format(output_file))

if __name__ == "__main__":
    main()
