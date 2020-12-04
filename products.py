#!/usr/bin/python3 -u

from classes.search_objects import ProductFinder
from classes.ucsl_io import read_blocks, find_ncs, make_block_stats, add_block
from classes.constants import Constants, ELB
import argparse
import os

script_path = os.path.split(os.path.realpath(__file__))[0]

def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("ncs", help='Specify NCS name or filename', type=str)
    parser.add_argument("elb_file", help='Specify ELB file', type=str)
    parser.add_argument("samples", help='Number of samples', type=int)
    parser.add_argument("--min_patterns", help='Minimal number of patterns', type=int, default=1)
    parser.add_argument("--output","-o", help='Output products to the file', type=str, default = '')
    return parser.parse_args()

def find_products(parameters, samples):
    scheme_optimizer = SchemeOptimizer(parameters)
    scheme_optimizer.run_products(samples)

def main():
    args = read_args()
    ncs, msg = find_ncs(args.ncs, script_path)
    if not ncs:
        print(msg)
        exit()
    result, blocks, ncs_name, deuterated = read_blocks(args.elb_file)
    if result:
        exit()
    product_finder = ProductFinder(blocks, args.samples, args.min_patterns)
    products = product_finder.find_products()
    schemes_total = sum([len(product) for product in products])
    blocks_total = {}
    print("{} products found".format(schemes_total))
    output = "[NCS = {}]\n".format(ncs_name)
    output+= "[Deuterated = {}]\n".format(deuterated)
    for product in products:
        for scheme in product:
            scheme.simplify()
            blocks_total = add_block(blocks_total, scheme.samples, len(scheme.patterns))
            output+=scheme.full_str()
            output+="\n"
    stats_total = make_block_stats(blocks_total)

    if args.output:
        with open(args.output, mode="w") as f:
            f.write(output)
            print("{} blocks were written to the file \'{}\'".format(schemes_total,args.output))
    else:
        print(output)
    print("Statistics of generated blocks:")
    print(stats_total["str"])



if __name__ == "__main__":
    main()




