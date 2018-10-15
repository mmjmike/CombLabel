#!/usr/bin/python3

from classes.interactive import *
from classes.constants import NCS, Constants, auto_name_ncs


DESCRIPTION = ("This script is designed to create NMR Coding System (NCS) file\n"
               "suitable for UUCSL package. To create a file you will need\n"
               "to specify amino acid isotopic labeling types (15N, 13C, 2H),\n"
               "NMR spectra that will be used and the name of the generated NCS\n"
               )


def input_deuteration():
    deuterated = False
    result = True
    query_text = "NMR coding system (NCS) will use DEUTERATED amino acids (Y[es]/N[o], default in No) > "
    while result:
        input_key = input(query_text)
        if answer_yes(input_key):
            deuterated = True
            print("Your are creating DEUTERATED NCS")
            break
        elif answer_no(input_key) or answer_empty(input_key):
            print("Your are creating NOT deuterated NCS")
            break
        else:
            result = ask_to_continue_input()
    return deuterated, result


def input_labels(deuterated):
    text = ("-------\n"
            "Please specify amino acid LABELING TYPES with one letter codes\n"
            "LABELING TYPE  X C A F N S T D  X F N D  <<< ONE-LETTER CODES\n"
            "    Alpha-15N  0 0 0 0 1 1 1 1  0 0 1 1\n"
            "        13-CO  0 1 0 1 0 1 0 1  0 1 0 1\n"
            "        13-CA  0 0 1 1 0 0 1 1  0 1 0 1\n"
            "   Deuterated  0 0 0 0 0 0 0 0  1 1 1 1\n"
            "\n"
            "Example: X N C D\n"
            )
    repeat_input = True
    result = False
    label_list = []
    while repeat_input:
        print(text)
        input_line = input("> ")
        result, label_list = strip_labels(input_line, deuterated)
        if result:
            print("you have specified following label types: ", ", ".join([label.name for label in label_list]))
            break
        else:
            # print("Please retry input")
            repeat_input = ask_to_continue_input()
    result = result and repeat_input
    return label_list, result


def strip_labels(line, deuterated):
    msg = ""
    label_types = []
    wrong_char = False
    for char in line:
        wrong_char = False
        if char.isalpha():
            if deuterated:
                if char.upper() not in Constants.DEUTERATED_TYPES:
                    msg = "Label type '{}' is not among deuterated label types".format(char)
                    wrong_char = True
            else:
                if char.upper() not in Constants.TYPES:
                    msg = "Character '{}' is not among possible label types".format(char)
                    wrong_char = True
            if wrong_char:
                print(msg)
                break
            else:
                label_types.append(Constants.BASIC_TYPES[Constants.TYPES_NAMES.index(char.upper())])
    return not wrong_char, label_types


def input_spectra():
    text = ("------\n"
            "Please specify spectra types with digits (1-{}), corresponding\n"
            "to the spectra, separated by spaces\n"
            "Warning: not including HSQC is not recommended\n".format(len(Constants.basic_spectra))
            )
    for i in range(len(Constants.basic_spectra)):
        spec = Constants.basic_spectra[i]
        text += "{}. {}\n".format(i+1, spec.name)
    text += "Example: 1 2 3"
    repeat_input = True
    result = False
    spectra_list = []
    while repeat_input:
        print(text)
        input_line = input("> ")
        result, spectra_list = strip_spectra(input_line)
        if result:
            print("you have specified following spectra: ", ", ".join([spectrum.name for spectrum in spectra_list]))
            break
        else:
            #print("Please retry input")
            repeat_input = ask_to_continue_input()
    result = result and repeat_input
    return spectra_list, result


def strip_spectra(line):
    number_list = line.split()
    num_list = set()
    spectra_list = []
    result = True
    for number in number_list:
        try:
            num = int(number)
        except ValueError:
            print("'{}' is not a digit".format(number))
            result = False
            break
        if num < 1 or num > len(Constants.basic_spectra):
            print("'{}' is not in range 1-{}".format(num,len(Constants.basic_spectra)))
            result = False
            break
        num_list.add(num-1)
    for i in range(len(Constants.basic_spectra)):
        spectrum = Constants.basic_spectra[i]
        if i in num_list:
            spectra_list.append(spectrum)
    return result, spectra_list


def input_name(ncs):
    repeat_input = True
    name = ''
    while repeat_input:
        print("===================================================")
        print("  This is the PREVIEW of your NCS:")
        print("===================================================")
        print(ncs)
        print("================= END OF PREVIEW ==================")
        auto_name = auto_name_ncs(ncs)
        text = "Please enter the NAME of the NCS above\n"
        text+= "Only letters, digits or symbols '-' and '_' can be used\n"
        text+= "We suggest the name \"{}\" for your NCS (Press enter to agree)".format(auto_name)
        input_line = input(text)
        good_name = True
        if answer_empty(input_line):
           good_name = True
           name = auto_name
           break
        for char in input_line:
            if not (char.isalpha() or char.isdigit() or char == "_" or char == "-"):
                print("Character '{}' can't be used in NCS name (not a digit or a letter)".format(char))
                good_name = False
                break
        if good_name:
            name = input_line
            print("Entered name is '{}'".format(name))
            break
        else:
            # print("Please repeat the input")
            repeat_input = ask_to_continue_input()
    return name, name


def print_description():
    print(DESCRIPTION)
    input("Press Enter to continue...")


def write_ncs_to_file(ncs):
    while True:
        input_query = "Do you want to save this NCS to the file \"{}.ncs\"? (Y[es]/N[o]) > ".format(ncs.name)
        input_key = input(input_query)
        if answer_yes(input_key) or answer_empty(input_key):
            with open("{}.ncs".format(ncs.name), "w") as f:
                f.write(str(ncs))
                print("The '{}' NCS was successfully written to '{}.ncs'".format(ncs.name, ncs.name))
                break
        elif answer_no(input_key):
            print("The NCS is NOT written to file")
            print("You can copy NCS to the file manually from the terminal output above")
            print("Bye")
            exit()


def main():
    print_description()
    deuterated_flag, continue_input = input_deuteration()
    if not continue_input:
        return
    labels, continue_input = input_labels(deuterated_flag)
    if not continue_input:
        return
    spectra, continue_input = input_spectra()
    if not continue_input:
        return
    ncs = NCS("NAME_OF_YOUR_NCS", spectra, labels, deuterated = deuterated_flag)
    name, continue_input = input_name(ncs)
    ncs.name = name
    if not continue_input:
        return
    write_ncs_to_file(ncs)


if __name__ == "__main__":
    main()
