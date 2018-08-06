#!/usr/bin/python3

from classes.constants import NCS, Constants

def answer_yes(input_key):
    return ( input_key.upper() == "Y" or input_key.upper() == "YES" )

def answer_no(input_key):
    return ( input_key.upper() == "N" or input_key.upper() == "NO" )

def answer_empty(input_key):
    return ( len(input_key) == 0  )


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


def ask_to_continue_input():
    result = False
    while True:
        input_key = input("Your answer is wrong. Do you want to try again? (Y[es]/n[o])> ")
        if answer_yes(input_key) or answer_empty(input_key):
            result = True
            break
        elif answer_no(input_key):
            print("See you later")
            exit()
    return result


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
            "Please specify spectra types with digits (1-7), corresponding\n"
            "to the spectra, separated by spaces\n"
            "Warning: not including HSQC is not recommended\n"
            "1. HSQC\n"
            "2. HNCO\n"
            "3. HNCA\n"
            "4. HNCOCA\n"
            "5. DQHNCA\n"
            "6. COHNCA\n"
            "7. HNCACO\n"
            "Example: 1 2 3"
            )
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
        if num < 1 or num > 7:
            print("'{}' is not in range 1-7".format(num))
            result = False
            break
        num_list.add(num-1)
    for i in range(len(Constants.basic_spectra)):
        spectrum = Constants.basic_spectra[i]
        if i in num_list:
            spectra_list.append(spectrum)
    return result, spectra_list

def auto_name_ncs(ncs):
    auto_name= ''
    if ncs.deuterated:
       auto_name+= "2H-"
    flag_X = False
    for type in ncs.label_types:
       if type.name == "X":
          flag_X = True
       if not ncs.deuterated and type.name == "X":
          continue
       auto_name+= type.name
    auto_name+= str(ncs.max_code_value)
    if not ncs.deuterated and not flag_X:
       auto_name+= "noX"
    return auto_name

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
    description = ("This script is designed to create NMR Coding System (NCS) file\n"
                   "suitable for UUCSL package. To create a file you will need\n"
                   "to specify amino acid isotopic labeling types (15N, 13C, 2H),\n"
                   "NMR spectra that will be used and the name of the generated NCS\n"
                   )
    print(description)
    input("Press Enter to continue...")


def write_to_file(ncs):
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
        
    print("===================================================")
    print("  This is again the PREVIEW of your NCS:")
    print("===================================================")
    print(ncs)
    print("================= END OF PREVIEW ==================")

    write_to_file(ncs)


if __name__ == "__main__":
    main()
