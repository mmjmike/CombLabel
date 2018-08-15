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
    parser.add_argument("file", help='Specify file for which to calculate stats', type=str)
    args = parser.parse_args()
    blocks = {}
    if not os.path.isfile(args.file):
        print("Error! File '{}' not found!".format(args.file))
        return
    print("File: '{}'".format(args.file))
    with open(args.file, "r") as f:
        for line in f:
            result = Constants.elb_re.match(line)
            if result:
                blocks = add_block(blocks, result.group(1), result.group(2))
    print(make_block_stats(blocks))


if __name__ == "__main__":
    main()