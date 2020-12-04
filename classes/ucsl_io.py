from cl_errors import errors as err
from classes.constants import Constants, NCS, ELB
# from classes.search_objects import ELB
import re
import os
import logging
import datetime
import random



def make_block_stats(blocks):
    output = "Blocks used:\n"
    total_blocks = 0
    samples = list(blocks.keys())
    samples.sort()
    blocks_by_samples = {}
    for samples_n in samples:
        pattern_numbers = list(blocks[samples_n].keys())
        pattern_numbers.sort()
        blocks_in_samples = 0
        for patterns_n in pattern_numbers:
            blocks_n = len(blocks[samples_n][patterns_n])
            total_blocks += blocks_n
            blocks_in_samples += blocks_n
            output += "{:2} X  {:2} - {} block(s)\n".format(samples_n, patterns_n, blocks_n)
        blocks_by_samples.update({samples_n: blocks_in_samples})
    output += "-----------\n"
    for samples_n in samples:
        output += "{:2} X  . - {} block(s)\n".format(samples_n, blocks_by_samples[samples_n])
    output += " TOTAL - {} blocks\n".format(total_blocks)
    output += "-----------\n"
    stats = {
        "str": output,
        "total": total_blocks
    }
    return stats


def read_blocks(block_file, logger=None, initial_blocks = {}):
    lines = read_lines(block_file)
    i = 0
    deuterated = False
    deuterated_found = False
    result = ''
    blocks = initial_blocks
    blocks_num = 0
    NCS_found = False
    ncs_name = ''
    ncs_regular = Constants.ncs_re
    while i < len(lines):
        if not NCS_found:
            ncs_match = ncs_regular.match(lines[i])
            if ncs_match:
                ncs_name = ncs_match.group(1).upper()
                NCS_found = True
        elif not deuterated_found:
            deuterated_match = Constants.deuterated_re.match(lines[i])
            if deuterated_match:
                deuteration = deuterated_match.group(1).upper()
                deuterated_found, deuterated = extract_deuterated(deuteration)
        else:
            elb_match = Constants.elb_re.match(lines[i])
            if elb_match:
                samples_num = int(elb_match.group(1))
                patterns_num = int(elb_match.group(2))
                blocks_num += 1
                sv_vec = []
                sv_match = Constants.sv_re.match(lines[i+1])
                if(sv_match):
                   sv_str = sv_match.group(1)
                   sv_vec = sv_str.split()
                   i = i + 1
                patterns = lines[i + 1:i + 1 + patterns_num]
                good_block = True
                for pattern in patterns:
                    if not good_block:
                        break
                    if len(pattern) != samples_num:
                        result = "Warning! The number of samples in block {} doesn't match " \
                              "\nwith specified value in blocks file '{}'.".format(blocks_num, block_file)
                        print(result)
                        if logger:
                            logger.warning(result)
                        good_block = False
                        break
                    for label_type in pattern:
                        if label_type not in Constants.TYPES:
                            result = "Warning! Unknown labeling type '{}' used in block {}" \
                                  "\n in blocks file '{}'.".format(label_type, blocks_num, block_file)
                            if logger:
                                logger.warning(result)
                            else:
                                print(result)
                            good_block = False
                            break
                if good_block:
                    block = ELB(patterns, ncs_name, deuterated, sv_vec)
                    if samples_num not in blocks:
                        blocks.update({samples_num: {patterns_num: [block]}})
                    else:
                        if patterns_num not in blocks[samples_num]:
                            blocks[samples_num].update({patterns_num: [block]})
                        else:
                            blocks[samples_num][patterns_num].append(block)
                i += patterns_num
        i += 1
    return result, blocks, ncs_name, deuterated


def add_block(blocks, samples, patterns):
    if samples not in blocks:
        blocks.update({samples: {patterns: [True]}})
    else:
        if patterns not in blocks[samples]:
            blocks[samples].update({patterns: [True]})
        else:
            blocks[samples][patterns].append(True)
    return blocks


def read_lines(filename):
    try:
        with open(filename, 'r', encoding='UTF-8') as f:
            lines = f.readlines()
    except IOError:
        raise err.ReadLinesError(filename)
    new_lines = []
    for line in lines:
        curr_line = line.rstrip()
        if curr_line != "" and curr_line[0] != "#":
            new_lines.append(curr_line.split("#")[0])
    return new_lines


def extract_labels(line):
    extracted_labels = [label.rstrip().upper() for label in line.split(",")]
    for label in extracted_labels:
        if label not in Constants.TYPES_NAMES:
            labels_names = ", ".join(Constants.TYPES_NAMES)
            print("Error! Label type '{}' is not in the following list:\n{}".format(label, labels_names))
            return []
    labels = []
    for label in Constants.BASIC_TYPES:
        if label.name in extracted_labels:
            labels.append(label)
    return labels


def extract_spectra(line):
    extracted_spectra = [spectrum.rstrip() for spectrum in line.split(",")]
    for spectrum in extracted_spectra:
        if spectrum not in Constants.SPECTRA_NAMES:
            spectra_names = ", ".join(Constants.SPECTRA_NAMES)
            print("Error! Spectrum '{}' is not in the following list:\n{}".format(spectrum, spectra_names))
            return []
    spectra = []
    for spectrum in Constants.basic_spectra:
        if spectrum.name in extracted_spectra:
            spectra.append(spectrum)
    return spectra


def extract_deuterated(deuteration):
    result = False
    deuterated = False
    d = deuteration.upper()
    if d == "TRUE" or d == "YES" or d == "Y" or d == "1":
        result = True
        deuterated = True
    elif d == "FALSE" or d == "NO" or d == "N" or d == "0":
        result = True
    return result, deuterated


def add_dir(path_list, new_dir):
    new_path_list = []
    for path in path_list:
        new_path_list.append(path)
        new_path_list.append(os.path.join(path, new_dir))
    return new_path_list


def add_file(path_list, filename):
    ncs_ext = ".ncs"
    new_path_list = []
    for path in path_list:
        new_path_list.append(os.path.join(path, filename))
        new_path_list.append(os.path.join(path, filename + ncs_ext))
    return new_path_list


def check_ncs_path(full_path):
    # print("Checking path: '{}'".format(full_path))
    if os.path.isfile(full_path):
        # print("File exists!")
        return True
    else:
        # print("File not found :(")
        return False


def add_to_file_name(filename, addition):
    path, only_name = os.path.split(filename)
    name, ext = os.path.splitext(only_name)
    name += addition
    only_name = name + ext
    return os.path.join(path, only_name)


def find_good_ncs_paths(ncs_name, script_path=os.path.split(os.path.realpath(__file__))[0]):
    ncs_dir = "NCS"
    paths = [os.getcwd(), script_path]
    paths = add_file(add_dir(paths, ncs_dir), ncs_name)
    good_paths = []
    for path in paths:
        if check_ncs_path(path):
            good_paths.append(path)
    return good_paths


def write_blocks(blocks, ncs_name, filename, deuterated, simplified = False):
    output = "[NCS = {}]\n".format(ncs_name)
    output+= "[Deuterated = {}]\n".format(deuterated)
    blocks_samples = list(blocks.keys())
    blocks_samples.sort()
    for samples_num in blocks_samples:
        patterns_numbers = list(blocks[samples_num].keys())
        patterns_numbers.sort()
        for patterns_num in patterns_numbers:
            for block in sorted(blocks[samples_num][patterns_num]):
                if(simplified):
                   block_head = "SIMPLIFIED"
                   block.simplify(join_char = "\n")
                   block_str = block.simplified_str
                else:
                   block_head = "ELB"
                   block_str = str(block)
                #output += "[{} samples = {} patterns = {}]\n".format(block_head, block.samples, len(block.patterns)) \
                #          + block_str
                output += block.full_str()+"\n"
    with open(filename, mode="w") as f:
        f.write(output)


def find_ncs(ncs_name, script_path):
    paths = find_good_ncs_paths(ncs_name, script_path)
    msg = ""
    for path in paths:
        msg += "Seek for NCS {} in path '{}'\n".format(ncs_name, path)
        ncs = read_ncs_file(path)
        if ncs:
            msg += "NCS '{}' read from '{}'\n".format(ncs_name, path)
            return ncs, msg
        else:
            msg += "NCS {} not found in path '{}'\n".format(ncs_name, path)
    msg += "Error! NCS file not found\n"
    return None, msg


def read_ncs_file(filename):
    name = ''
    deuterated = False
    deuterated_read = False
    labels = []
    spectra = []
    with open(filename, "r") as f:
        for line in f:
            if not name:
                name_match = Constants.ncs_re.match(line)
                if name_match:
                    name = name_match.group(1)
                    continue
            if not labels:
                labels_match = Constants.labels_re.match(line)
                if labels_match:
                    labels = extract_labels(labels_match.group(1))
                    continue
            if not spectra:
                spectra_match = Constants.spectra_re.match(line)
                if spectra_match:
                    spectra = extract_spectra(spectra_match.group(1))
                    continue
            if not deuterated_read:
                deuterated_match = Constants.ncs_re.match(line)
                if deuterated_match:
                    deuteration = deuterated_match.group(1).upper()
                    deuterated_found, deuterated = extract_deuterated(deuteration)
                    continue
    if name and labels and spectra:
        return NCS(name, spectra, labels, deuterated)
    return None


def read_sequence(filename):
    sequence = ""
    lines = read_lines(filename)
    for line in lines:
        if line[0] == ">":
            continue
        clean_line = line.replace("\t", "").replace(" ", "")
        for char in line:
            if char.upper() not in Constants.RES_TYPES_LIST:
                print("Error! Symbol '{}' doesn't represent any amino acid residue in sequence file '{}'. Please correct sequence file".format(char, filename))
                return ""
        sequence += clean_line.upper()
    return sequence


def read_solution(filename):
    solution_dict = {}
    lines = read_lines(filename)
    solution_start = False
    header_read = False
    samples = 0
    for line in lines:
        solution_match = Constants.solution_re.match(line)
        if solution_match:
            solution_start = True
            continue
        if solution_start:
            if not header_read:
                header = [value.strip() for value in line.split(",")]
                if header[0].upper() == "RES" and len(header) > 2:
                    header_read = True
                    samples = len(header) - 1
                    continue
            else:
                data = [value.strip() for value in line.split(",")]
                if data[0].capitalize() not in Constants.RES_TYPES_THREE:
                    break
                if len(data) != samples + 1:
                    print("Error in solution file '{}'! Length of line '{}' is not the same as the header length".format(filename, line))
                    return {}
                if data[0].upper() in Constants.RES_TYPES_LIST:
                    residue = data[0].upper()
                elif data[0].capitalize() in Constants.RES_TYPES_THREE:
                    residue = Constants.TO_ONE_LETTER_CODE[data[0].capitalize()]
                else:
                    print(
                        "Error in solution file '{}'! Unknown residue type '{}'".format(
                            filename, data[0]))
                    return {}
                pattern = ""
                for label in data[1:]:
                    if label.upper() in Constants.TYPES:
                        pattern += label.upper()
                    else:
                        print(
                            "Error in solution file '{}'! Unknown label type '{}'".format(
                                filename, label))
                        return {}
                solution_dict[residue] = pattern
    return solution_dict


def write_ncs_stamp(ncs):
    return "[NCS = {}]\n[Deuterated = {}]".format(ncs.name, ncs.deuterated)


def write_best_scheme(best_scheme, filename):
    output = '________________________\nBest scheme\n[NCS = {}]\n'.format(best_scheme.scheme.ncs_name)
    output += "[solution]\nRes"
    for i in range(best_scheme.samples):
        output += ",S{}".format(i + 1)
    output += "\n"
    for res in best_scheme.residues:
        output += res + ", " + ", ".join(list(best_scheme.label_dict[res])) + "\n"
    output += "[price]\n{}\n\n".format(best_scheme.price)
    output += "Blocks used for this scheme:\n"
    for block in best_scheme.blocks:
        output += block.full_str() + "\n"
    with open(filename, mode="w") as f:
        f.write(output)
    return output


def write_product_stats(stats, filename):
    output = ""
    for samples in sorted(stats):
        output += "-------------\n"
        output += "Products in {} samples:\n".format(samples)
        for product_type in sorted(stats[samples]):
            output += product_type + "\n"
            for status in sorted(stats[samples][product_type]):
                output += "{:4} scheme(s): {}\n".format(stats[samples][product_type][status], status)
    with open(filename, mode="w") as f:
        f.write(output)
    return output


def write_products(products, samples, filename, mode='w'):
    if not products:
        return
    product_schemes = 0
    output = ""
    # if mode == 'w':
    #     output += write_ncs_stamp(products[0].ncs)
    output += make_block_stats(products[0].blocks)["str"]
    output += "\n-----------------------\nProducts calculated for {} samples:\n".format(samples)
    for product in products:
        output += str(product) + ": {} scheme(s)\n".format(len(product))
        product_schemes += len(product)
    output += "{} product types calculated\n".format(len(products))
    output += "{} total labeling schemes to check\n".format(product_schemes)
    with open(filename, mode=mode) as f:
        f.write(output)
    return output


def read_prices(prices_file):
    msg = ""
    deuterated_found = False
    prices = {}
    lines = read_lines(prices_file)
    if len(lines) < 3:
        msg = "Prices file '{}' is too short".format(prices_file)
        return prices, msg
    d = []
    line_length = 0
    first_line = True
    for line in lines:
        deuterated_match = Constants.deuterated_re.match(line)
        if deuterated_match:
            deuteration = deuterated_match.group(1).upper()
            deuterated_found, deuterated = extract_deuterated(deuteration)
            prices.update({"Deuterated" : deuterated})
            continue
        else:
          s = [x.strip() for x in line.split(",")]
          if first_line:
              line_length = len(s)
              first_line = False
          else:
              if len(s) != line_length:
                  msg = "Not equal length of lines in prices file '{}'".format(prices_file)
                  return prices, msg
          d.append(s)
    price_label_types = d[0][1:]
    for l_type in price_label_types:
        if l_type.upper() not in Constants.TYPES_NAMES:
            msg = "Incorrect label type {} in prices file '{}'".format(l_type, prices_file)
            return prices, msg
    if len(price_label_types) < 2:
        msg = "Too few labeling types specified in prices file '{}'".format(prices_file)
        return prices, msg

    for i in range(len(d) - 1):
        curr_dict = {}
        for j in range(len(price_label_types)):
            try:
                price = float(d[i + 1][j + 1])
            except ValueError:
                msg = "ERROR! Price must be set in digits"
                msg += "\nPlease check price file '{}' (row {}; col {})".format(prices_file,
                                                                                i + 2, j + 2)
                return {}, msg
            curr_dict.update({price_label_types[j]: price})
        residue_type = d[i + 1][0]
        if residue_type not in Constants.RES_TYPES_LIST and residue_type not in Constants.RES_TYPES_THREE:
            msg = "Wrong residue type '{}' in prices file '{}'".format(residue_type, prices_file)
            return {}, msg
        if residue_type in Constants.RES_TYPES_THREE:
            residue_type1 = Constants.TO_ONE_LETTER_CODE[residue_type]
            residue_type = residue_type1
        prices.update({residue_type: curr_dict})
    return prices, msg
