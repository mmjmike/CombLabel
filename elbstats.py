#!/usr/bin/python3 -u

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
    parser.add_argument("--sv", "-s", help='Print SV min andd max values', action="store_true")
    args = parser.parse_args()
    blocks_total={}
    sv_min = dict()
    sv_max = dict()
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
                    samples = int(result.group('samples'))
                    patterns = int(result.group('patterns'))
                    #print("{:>2} {:>2}".format(samples, patterns))
                    blocks = add_block(blocks, samples, patterns)
                    blocks_total = add_block(blocks_total, samples, patterns)
                    continue
                if args.sv:
                    sv_re_match = Constants.sv_re.match(line)
                    if sv_re_match:
                       sv = [int(i) for i in sv_re_match.group('sv').split()]
                       if samples in sv_min.keys():
                           sv_min[samples] = list(map(min, zip(sv_min[samples], sv))) 
                           sv_max[samples] = list(map(max, zip(sv_max[samples], sv))) 
                       else:
                           #print(sv)
                           #print(sv_min.keys())
                           sv_min[samples] = sv
                           sv_max[samples] = sv
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
    if args.sv:
        print("Min and max values for SV vectors:")
        for s in sorted(sv_min.keys()):
            
            min_out = "SV_MIN ({} samples): ".format(s)
            max_out = "SV_MAX ({} samples): ".format(s)
            #print(sv_min)
            #print(sv_max)
            for mins in sv_min[s]:
                min_out += "{:>2} ".format(mins)
            for maxs in sv_max[s]:
                max_out += "{:>2} ".format(maxs)
            print(min_out)
            print(max_out)

if __name__ == "__main__":
    main()
