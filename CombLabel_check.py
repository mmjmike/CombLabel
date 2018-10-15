#!/usr/bin/python3

import argparse
import os
from classes.ucsl_io import read_sequence, read_solution, find_ncs, read_prices
from classes.constants import Constants, NCS, auto_name_ncs
from classes.search_objects import PairsTable, Sequence, calc_price



def get_stock_from_solution(solution):
    label_options = {}
    for residue in solution:
        labels = []
        pattern = solution[residue]
        for label_type in Constants.BASIC_TYPES:
            if label_type.name in pattern:
                labels.append(label_type)
        label_options[residue] = labels
    return label_options


def add_unlabeled_residues(solution, sequence):
    res_types = sequence.residue_types_used
    samples = len(next(iter(solution.values())))
    for residue in res_types:
        if residue not in solution:
            solution[residue] = "X" * samples
    return solution


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

    protein_sequence = read_sequence(args.sequence)
    if len(protein_sequence) == 0:
        exit()
    if len(protein_sequence) < 3:
        print("Protein sequence is too short({})".format(len(protein_sequence)))
        exit()
    sequence = Sequence(protein_sequence)
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
    output_solution = name + "_solution.txt"
    output_dictionary = name + "_dictionary.txt"

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


def write_solution_stats(parameters):
    codes, meanings, meanings_dict = calculate_code_dictionary(parameters)
    write_code_dict(parameters, codes, meanings)

    stats = calculate_stats(parameters, meanings_dict)
    pairs_table = parameters["pairs_table"]

    output = ""
    output += pairs_table.make_full_pairs_table()
    output += pairs_table.make_stock_pairs_table()
    output += str(parameters["ncs"])
    output += make_solution_output(parameters)
    output += pairs_table.make_pairs_codes(parameters["solution"], parameters["ncs"])
    output += stats_to_text(stats)
    with open(parameters["output_solution"], "w") as f:
        f.write(output)
        f.close()
    print("Solution stats written to file '{}'".format(parameters["output_solution"]))


def calculate_code_dictionary(parameters):
    codes_table = []
    meanings_table = []
    code_meanings = {}
    sequence = parameters["sequence"]
    solution = parameters["solution"]
    ncs = parameters["ncs"]

    for i in range(len(sequence.sequence) - 1):
        pair = sequence.sequence[i:i+2]
        residue1 = pair[0]
        residue2 = pair[1]
        pattern1 = solution[residue1]
        pattern2 = solution[residue2]
        code = ncs.calc_code(pattern1, pattern2)
        meaning = "".join([residue1, str(i+1), " - ", residue2, str(i+2)])
        codes_table.append(code)
        meanings_table.append(meaning)
        if code not in code_meanings:
            code_meanings[code] = [pair]
        else:
            code_meanings[code].append(pair)

    sorted_meanings = [mean for code, mean in sorted(zip(codes_table,
                                                            meanings_table))]
    sorted_codes_table = [code for code in sorted(codes_table)]
    return sorted_codes_table, sorted_meanings, code_meanings


def write_code_dict(parameters, codes, meanings):
    output = "#"*50+"\n"
    output += "# The Spectrum codes dictionary for \'"+parameters["name"]+"\':\n"
    output += "# Spectrum code: First AA - Second AA\n\n"
    for i in range(len(codes)):
        output += "{}: {}\n".format(codes[i], meanings[i])
    with open(parameters["output_dictionary"], "w") as f:
        f.write(output)
        f.close()
    print("Code dictionary written to file '{}'".format(parameters["output_dictionary"]))


def make_solution_output(parameters):
    solution = parameters["solution"]
    pairs_table = parameters["pairs_table"]
    output = "[solution]\n"
    output += "% Soluiton number = 0\n"
    if parameters["show_price"]:
        price = 0
        for res in solution:
            price += calc_price(parameters["price_dict"], res, solution[res])
        output += "% Solution price  = {}".format(str(round(price, 2)))
    output += "Res"
    samples_num = len(next(iter(solution.values())))
    for i in range(samples_num):
        output += ",S" + str(i + 1)
    output += "\n"
    for residue in pairs_table.residues_to_label:
        output += Constants.TO_THREE_LETTER_CODE[residue] + ", " + ", ".join(
            list(solution[residue])) + "\n"
    for res in pairs_table.non_labeled_residues:
        output += Constants.TO_THREE_LETTER_CODE[res] + ", " + ", ".join(list("X" * samples_num)) + "\n"
    output += "\n"
    return output


def calculate_stats(parameters, meanings_dict):
    stats = {}
    sequence = parameters["sequence"]
    pairs_table = parameters["pairs_table"]

    N_aa = len(sequence.sequence)
    N_pairs = len(pairs_table.unique_pairs)
    PI_p = len([pair for pair in pairs_table.unique_pairs if pair[1] == "P"])
    PI_a = len([res for res in sequence.sequence if res == 'P'])
    PU_p = 0
    PN_a = 0
    PN_p = 0
    for i in range(len(pairs_table.residues_first)):
        for j in range(len(pairs_table.residues_second)):
            number_of_pairs = pairs_table.all_residue_pairs[i][j]
            if number_of_pairs == 1:
                PU_p += 1
            elif number_of_pairs > 1:
                PN_p += 1
                PN_a += number_of_pairs
    PU_a = PU_p
    SI_a = 0
    SI_p = 0
    SU_p = 0
    SN2_a = 0
    SN2_p = 0
    SN1_a = 0
    SN1_p = 0
    for i in range(len(pairs_table.residues_carbon)):
        for j in range(len(pairs_table.residues_nitro)):
            cell_value = pairs_table.residue_pairs[i][j]
            res1 = pairs_table.residues_carbon[i]
            res2 = pairs_table.residues_nitro[j]
            if cell_value:
                seq = sequence.sequence
                if res2 == 'Other':
                    diff_pairs = len(set([seq[s + 1] for s in range(len(seq) - 1)
                                          if seq[s] == res1
                                          and seq[s + 1] in pairs_table.residues_not_nitro]))
                    if res1 == 'Other':
                        diff_pairs = len(set([seq[s + 1] for s in range(len(seq) - 1)
                                              if seq[s] in pairs_table.residues_not_carbon
                                              and seq[s + 1] in pairs_table.residues_not_nitro]))
                    SI_a += cell_value
                    SI_p += diff_pairs
                elif res1 == 'Other':
                    diff_pairs = len(set([seq[s] for s in range(len(seq) - 1)
                                          if seq[s + 1] == res2
                                          and seq[s] in pairs_table.residues_not_carbon]))
                    if res2 in pairs_table.bad_residues:
                        SN1_p += diff_pairs
                        SN1_a += cell_value
                    elif cell_value == 1:
                        SU_p += 1
                    elif diff_pairs > 1:
                        SN1_p += diff_pairs
                        SN1_a += cell_value
                    else:
                        SN2_p += 1
                        SN2_a += cell_value
                elif (res1 == res2
                      and res2 in pairs_table.bad_residues):
                    SN1_p += 1
                    SN1_a += cell_value
                elif cell_value == 1:
                    SU_p += 1
                else:
                    SN2_p += 1
                    SN2_a += cell_value
    SU_a = SU_p
    stats.update({
        "SIa": SI_a,
        "SIp": SI_p,
        "SUa": SU_a,
        "SUp": SU_p,
        "SN1a": SN1_a,
        "SN1p": SN1_p,
        "SN2a": SN2_a,
        "SN2p": SN2_p,
        "Na": N_aa,
        "Np": N_pairs,
        "PIa": PI_a,
        "PIp": PI_p,
        "PUa": PU_a,
        "PUp": PU_p,
        "PNa": PN_a,
        "PNp": PN_p
    })

    LI_a = 0
    LI_p = 0
    LN1_a = 0
    LN1_p = 0
    LN2_a = 0
    LN2_p = 0
    LU_a = 0
    LA_p = 0
    LA_a = 0

    for code in meanings_dict:
        meaning = meanings_dict[code]
        if set(list(code)) == set("0"):
            LI_p += len(set(meaning))
            LI_a += len(meaning)
        elif len(meaning) == 1:
            LU_a += 1
        elif len(set(meaning)) == 1:
            LN2_p += 1
            LN2_a += len(meaning)
        elif len(set([pair[1] for pair in meaning])) == 1:  # all second residues in pair are equal
            LN1_p += len(set(meaning))
            LN1_a += len(meaning)
        else:
            LA_a += len(meaning)
            LA_p += len(set(meaning))
    LU_p = LU_a

    stats.update({
        "LIa": LI_a,
        "LIp": LI_p,
        "LUa": LU_a,
        "LUp": LU_p,
        "LN1a": LN1_a,
        "LN1p": LN1_p,
        "LN2a": LN2_a,
        "LN2p": LN2_p,
        "LAa": LA_a,
        "LAp": LA_p
    })
    return stats


def stats_to_text(stats):
    output = "\n\n" + "#" * 50 + "\n"
    output += "# Calculation statistics\n"
    output += "# \n\n"
    output += "[stats]\n"
    output += ""
    output += "\n\n"
    output += "# Statistics for PAIRS in amino acid sequence\n"
    output += "# The STOCK (availability of isotopically labeled amino acid)\n"
    output += "# is NOT accounted for in this statistics\n"
    output += "# The labeling scheme is NOT accounted too\n"
    output += "\n"
    output += "[stats,pairs]\n"
    output += "Par, Description,  Residues, Pairs\n"
    output += "N,   Total number, {:>8}, {:>5}\n".format(stats["Na"], stats["Np"])
    output += "PI,  Invisible,    {:>8}, {:>5}\n".format(stats["PIa"], stats["PIp"])
    output += "PU,  Unique,       {:>8}, {:>5}\n".format(stats["PUa"], stats["PUp"])
    output += "PN,  Non-unique,   {:>8}, {:>5}\n".format(stats["PNa"], stats["PNp"])
    output += "\n"
    output += "# Statistics for STOCK-available pairs in amino acid sequence\n"
    output += "# The STOCK is used to check whether the particular pairs are distinguishable \n"
    output += "# in principle with some labeling scheme unlimited in size with some NMR spectra\n"
    output += "# The particular labeling scheme, found by the program, is NOT accounted here\n"
    output += "\n"
    output += "[stats,stock]\n"
    output += "Par, Description,     Residues, Pairs\n"
    output += "N,   Total number,    {:>8}, {:>5}\n".format(stats["Na"], stats["Np"])
    output += "SI,  Invisible,       {:>8}, {:>5}\n".format(stats["SIa"], stats["SIp"])
    output += "SU,  Unique code,     {:>8}, {:>5}\n".format(stats["SUa"], stats["SUp"])
    output += "SN2, AA type of both, {:>8}, {:>5}\n".format(stats["SN2a"], stats["SN2p"])
    output += "SN1, AA type of last, {:>8}, {:>5}\n".format(stats["SN1a"], stats["SN1p"])
    output += "\n"
    output += "# Statistics for LABELING CODES\n"
    output += "# The pairs are distinguishable, if their labeling codes are different\n"
    output += "# Both sequence, stock, NMR spectra and particular labeling scheme is accounted here\n"
    output += "\n"
    output += "[stats,labeling]\n"
    output += "Par, Description,     Residues, Pairs\n"
    output += "N,   Total number,    {:>8}, {:>5}\n".format(stats["Na"], stats["Np"])
    output += "LI,  Invisible,       {:>8}, {:>5}\n".format(stats["LIa"], stats["LIp"])
    output += "LU,  Unique code,     {:>8}, {:>5}\n".format(stats["LUa"], stats["LUp"])
    output += "LN2, AA type of both, {:>8}, {:>5}\n".format(stats["LN2a"], stats["LN2p"])
    output += "LN1, AA type of last, {:>8}, {:>5}\n".format(stats["LN1a"], stats["LN1p"])
    output += "LA,  Ambiguous code,  {:>8}, {:>5}\n".format(stats["LAa"], stats["LAp"])
    return output


def main():

    args = get_args()
    parameters = get_params(args)
    write_solution_stats(parameters)


if __name__ == "__main__":
    main()