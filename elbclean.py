#!/usr/bin/python3 -u

from classes.search_objects import ProductFinder
from classes.ucsl_io import make_block_stats, add_to_file_name, write_blocks, read_blocks
import argparse
import os

# add logging, especially in the reading of blocks
LOG_ITERATION = 10000


def clear_identical_blocks(blocks):
    iter = 0
    total_blocks = make_block_stats(blocks)["total"]
    for samples_num in blocks:
        patterns_numbers = list(blocks[samples_num].keys())
        patterns_numbers.sort(reverse=True)
        for pattern_num in patterns_numbers:
            new_block_dict = {}
            for block_one in blocks[samples_num][pattern_num]:
                iter += 1
                if iter % LOG_ITERATION == 0:
                    print("Identical blocks checked {}/{}".format(iter, total_blocks))
                new_block_dict[block_one] = block_one
            if not new_block_dict:
                blocks[samples_num].pop(pattern_num)
            else:
                blocks[samples_num][pattern_num] = list(new_block_dict.values())
    return clear_empty_block_types(blocks)

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
                    print("Redundant blocks checked {}/{}".format(iter, total_blocks))
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
                    print("Product blocks checked {}/{}".format(iter, total_blocks))
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


def clear_blocks(blocks, flags):
    (identical_flag, redundant_flag, product_flag) = flags
    print(make_block_stats(blocks)["str"]+"\n")

    if identical_flag:
        print("Clearing identical blocks...")
        blocks = clear_identical_blocks(blocks)
        print("Identical blocks cleared...")
        print(make_block_stats(blocks)["str"]+"\n")

    if redundant_flag:
        print("Clearing redundant blocks...")
        blocks = clear_redundant_blocks(blocks)
        print("Redundant blocks cleared...")
        print(make_block_stats(blocks)["str"]+"\n")

    if product_flag:
        print("Clearing product blocks...")
        blocks = clear_product_blocks(blocks)
        print("Product blocks cleared...")
        print(make_block_stats(blocks)["str"]+"\n")
    return blocks


def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("elb_files", help='Specify ELB file(s)', type=str, nargs='+')
    parser.add_argument("-o", dest="output_file", help='Specify output file (optional)', default="", type=str)
    parser.add_argument("--identical", "-I", help='Clear identical blocks [Default]', action="store_true", default = True)
    parser.add_argument("--identical-skip", "-i", dest="identical", help='Do NOT clear identical blocks', action="store_false")
    parser.add_argument("--redundant", "-R", help='Clear redundant blocks [Default]', action="store_true", default = True)
    parser.add_argument("--redundant-skip", "-r", dest="redundant", help='Do NOT clear redundant blocks', action="store_false")
    parser.add_argument("--product", "-P", help='Clear product blocks [Default]', action="store_true", default = True)
    parser.add_argument("--product-skip", "-p", dest="product", help='Do NOT clear product blocks', action="store_false")
    parser.add_argument("--simplified", "-s", dest="simplified", help='Output SIMPLIFIED blocks', action="store_true")
    args = parser.parse_args()
    print("Clear IDENTICAL blocks: {}".format(args.identical))
    print("Clear REDUNDANT blocks: {}".format(args.redundant))
    print("Clear PRODUCT blocks:   {}".format(args.product))
    print("Output SIMPLIFIED blocks:   {}".format(args.simplified))
    for file in args.elb_files:
        if not os.path.isfile(file):
            print("Error! ELB file '{}' not found".format(file))
            exit()
        print("ELB file: {}".format(file))
    if args.output_file:
        print("Output file: {}".format(args.output_file))
        output_file = args.output_file
    else:
        output_file = add_to_file_name(args.elb_files[0] + "_" + args.elb_files[-1], "_clean")
    return args.elb_files, output_file, (args.identical, args.redundant, args.product), args.simplified


def main():
    elb_files, output_file, flags, output_simplified = read_args()
    blocks = {}
    for file in elb_files:
        result, blocks, ncs_name, deuterated = read_blocks(file, initial_blocks = blocks)
        print("{} read".format(file))
    if result:
        print(result)
        return
    new_blocks = clear_blocks(blocks, flags)
    write_blocks(new_blocks, ncs_name, output_file, deuterated, simplified = output_simplified )
    print("Blocks written to '{}' successfully".format(output_file))


if __name__ == "__main__":
    main()
