#!/usr/bin/python3

from classes.search_objects import ProductFinder
from classes.ucsl_io import make_block_stats, add_to_file_name, write_blocks, read_blocks
import argparse
import os
import logging

# add logging, especially in the reading of blocks
LOG_ITERATION = 10000


def clear_redundant_blocks(blocks):
    iter = 0
    total_blocks = make_block_stats(blocks)["total"]
    for samples_num in blocks:
        patterns_numbers = list(blocks[samples_num].keys())
        patterns_numbers.sort(reverse=True)
        good_blocks = []
        for pattern_num in patterns_numbers:
            new_block_list = []
            for block in blocks[samples_num][pattern_num]:
                iter += 1
                if iter % LOG_ITERATION == 0:
                    print("Blocks checked {}/{}".format(iter, total_blocks))
                block_good = True
                for good_block in good_blocks:
                    if block.is_subset_of(good_block.simplified):
                        block_good = False
                        break
                if block_good:
                    good_blocks.append(block)
                    new_block_list.append(block)
            if not new_block_list:
                blocks[samples_num].pop(pattern_num)
            else:
                blocks[samples_num][pattern_num] = new_block_list
    return clear_empty_block_types(blocks)


def clear_product_blocks(blocks):
    blocks_samples = list(blocks.keys())
    blocks_samples.sort(reverse=True)
    iter = 0
    total_blocks = make_block_stats(blocks)["total"]
    for samples_num in blocks_samples:
        patterns_numbers = list(blocks[samples_num].keys())
        patterns_numbers.sort(reverse=True)
        for patterns_num in patterns_numbers:
            good_blocks = []
            for block in blocks[samples_num][patterns_num]:
                iter += 1
                if iter % LOG_ITERATION == 0:
                    print("Blocks checked {}/{}".format(iter, total_blocks))
                block_good = True
                prod_finder = ProductFinder(blocks, samples_num, patterns_num, equal=True)
                products = prod_finder.find_products()
                for product in products:
                    if len(product.product_list) == 1:
                        continue
                    if not block_good:
                        break
                    for product_block in product:
                        if product_block == block:
                            block_good = False
                            break
                if block_good:
                    good_blocks.append(block)
            blocks[samples_num][patterns_num] = good_blocks
    return clear_empty_block_types(blocks)


def clear_empty_block_types(blocks):
    bad_samples_keys = []
    for samples_num in blocks:
        bad_patt_keys = []
        for patterns_num in blocks[samples_num]:
            if not blocks[samples_num][patterns_num]:
                bad_patt_keys.append(patterns_num)
        for patterns_num in bad_patt_keys:
            blocks[samples_num].pop(patterns_num)
        if blocks[samples_num] == {}:
            bad_samples_keys.append(samples_num)
    for samples_num in bad_samples_keys:
        blocks.pop(samples_num)
    return blocks


def clear_blocks(blocks):
    print(make_block_stats(blocks)["str"])
    print("Clearing redundant blocks...\n")
    blocks = clear_redundant_blocks(blocks)
    print("Redundant blocks cleared...\n")
    print(make_block_stats(blocks)["str"])
    print("Clearing product blocks...\n")
    blocks = clear_product_blocks(blocks)
    print("Product blocks cleared...\n")
    print(make_block_stats(blocks)["str"])
    return blocks


def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("elb_file", help='Specify ELB file', type=str)
    parser.add_argument("output_file",
                        help='Specify output file (optional)', default="",
                        type=str, nargs="?")
    args = parser.parse_args()
    if not os.path.isfile(args.elb_file):
        print("Error! ELB file '{}' not found".format(args.elb_file))
        exit()
    print("ELB file: {}".format(args.elb_file))
    if args.output_file:
        print("Output file: {}".format(args.output_file))
        output_file = args.output_file
    else:
        output_file = add_to_file_name(args.elb_file, "_clean")
    return args.elb_file, output_file


def main():
    elb_file, output_file = read_args()
    result, blocks, ncs_name, deuterated = read_blocks(elb_file)
    if result:
        print(result)
        return
    new_blocks = clear_blocks(blocks)
    write_blocks(new_blocks, ncs_name, output_file, deuterated)
    print("Blocks written to '{}' successfully".format(output_file))


if __name__ == "__main__":
    main()
