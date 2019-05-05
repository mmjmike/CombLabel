import argparse
import os
from classes.ucsl_io import read_sequence, read_solution, find_ncs, read_prices
from classes.constants import Constants, NCS, auto_name_ncs
from classes.search_objects import CLSequence
from classes.sequence_specific import PairsTable, write_solution_stats,\
    add_unlabeled_residues, get_stock_from_solution


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("ncs", help='Specify NCS (NMR Coding System)', type=str)
    parser.add_argument("sequence",
                        help='Specify sequence of protein of interest in one-letter amino acid code ()',
                        type=str)
    parser.add_argument("solution",
                        help='Specify solution file',
                        type=str)
    parser.add_argument('--prices', '-p',
                        help='Specify prices file for assessment of solution price', default="",
                        type=str)
    return parser.parse_args()


def get_params(args):
    script_path = os.path.split(os.path.realpath(__file__))[0]
    ncs, msg = find_ncs(args.ncs, script_path)
    if not ncs:
        print(msg)
        exit()
    print("NCS '{}' read".format(ncs.name))

    protein_sequence, msg = read_sequence(args.sequence)
    print("Protein sequence: {}".format(protein_sequence))
    if len(protein_sequence) == 0:
        exit()
    if len(protein_sequence) < 3:
        print("Protein sequence is too short({})".format(len(protein_sequence)))
        exit()
    sequence = CLSequence(protein_sequence)
    print("Protein sequence file '{}' read".format(args.sequence))

    solution = read_solution(args.solution)
    if not solution:
        exit()
    solution = add_unlabeled_residues(solution, sequence)
    print("Solution file '{}' read".format(args.solution))

    label_options = get_stock_from_solution(solution)
    pairs_table = PairsTable(sequence, label_options)
    all_label_options = set(sum(list(label_options.values()), []))

    new_label_types = []
    for label_type in Constants.BASIC_TYPES:
        if label_type in all_label_options or label_type in ncs.label_types:
            new_label_types.append(label_type)

    show_price = False
    price_dict = {}
    if args.prices:
        show_price = True
        price_dict, msg = read_prices(args.prices)
        if not price_dict:
            print(msg)
            exit()


    new_ncs = NCS("", ncs.spec_list, new_label_types, deuterated=ncs.deuterated)
    new_ncs.name = auto_name_ncs(new_ncs)


    name = "{}_{}".format(args.sequence.split(".")[0], new_ncs.name)
    output_solution = name + "_solution_c.txt"
    output_dictionary = name + "_dictionary_c.txt"

    params = {
        "ncs": new_ncs,
        "name": name,
        "label_options": label_options,
        "sequence": sequence,
        "solution": solution,
        "pairs_table": pairs_table,
        "show_price": show_price,
        "price_dict": price_dict,
        "output_solution": output_solution,
        "output_dictionary": output_dictionary
    }
    return params


def main():

    args = get_args()
    parameters = get_params(args)
    write_solution_stats(parameters)


if __name__ == "__main__":
    main()