#!/usr/bin/python3 -u

from classes.search_objects import ProductFinder
from classes.ucsl_io import make_block_stats, add_to_file_name, write_blocks, read_blocks
import argparse
import os

def pattern_stock(blocks):
    ps_list=[]
    for samples_num in blocks:
        pset=set()
        patterns_numbers = list(blocks[samples_num].keys())
        patterns_numbers.sort(reverse=True)
        for pattern_num in patterns_numbers:
            for block in blocks[samples_num][pattern_num]:
                for p in block.patterns:
                    pset.add(p)
        pl=list(pset)
        print(pl)
        pl.sort()
        print(pl)
        ps_list.append(pl)
    return ps_list

def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("elb_files", help='Specify ELB file(s)', type=str, nargs='+')
    args = parser.parse_args()
    print(args.elb_files)
    return args.elb_files

def main():
    elb_files = read_args()
    blocks = {}
    for file in elb_files:
        result, blocks, ncs_name, deuterated = read_blocks(file, initial_blocks = blocks)
        print("{} read".format(file))

    pslist = pattern_stock(blocks)
    print("pslist generated")
    for s in pslist:
        print(s)
        for p in s:
            print(p)

if __name__ == "__main__":
    main()