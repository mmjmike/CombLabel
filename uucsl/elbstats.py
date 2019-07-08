#!/usr/bin/python3

import argparse, os
from classes.ucsl_io import make_block_stats
from classes.constants import Constants


def add_block(blocks, samples, patterns):
    if samples not in blocks:
        blocks.update({samples: {patterns: [True]}})
    else:
        if patterns not in blocks[samples]:
            blocks[samples].update({patterns: [True]})
        else:
            blocks[samples][patterns].append(True)
    return blocks


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("files", help='Specify file(s) for which to calculate stats', type=str, nargs='+')
    parser.add_argument("--total", "-t", help='Print only total number of blocks', action="store_true")
    args = parser.parse_args()
    blocks_total={}
    for file in args.files:
        blocks = {}
        if not os.path.isfile(file):
            print("Error! File '{}' not found!".format(file))
            return
        output = "{:>25}:".format(file)
        with open(file, "r") as f:
            for line in f:
                result = Constants.elb_re.match(line)
                if result:
                    blocks = add_block(blocks, int(result.group(1)), int(result.group(2)))
                    blocks_total = add_block(blocks_total, int(result.group(1)), int(result.group(2)))
        stats = make_block_stats(blocks)
        if args.total:
            output += " {:>7}".format(int(stats["total"]))
        else:
            output += "\n" + stats["str"]
        print(output)
    stats_total = make_block_stats(blocks_total)
    if len(args.files)>1:
        output = "{:>25}:".format("TOTAL FOR {} FILES".format(len(args.files)))
        if args.total:
            output += " {:>7}".format(int(stats_total["total"]))
        else:
            output += "\n" + stats_total["str"]
        print(output)

if __name__ == "__main__":
    main()
