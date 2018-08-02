from classes.constants import NCS, Constants


def input_deuteration():
    deuterated = False
    result = True
    while result:
        input_key = input("Are the amino acid label types deuterated? (Y/N) > ")
        if input_key.upper() == "Y":
            deuterated = True
            break
        elif input_key.upper() == "N":
            break
        else:
            result = ask_to_continue_input()
    return deuterated, result


def ask_to_continue_input():
    result = False
    while True:
        input_key = input("Want to continue? (Y/N)> ")
        if input_key.upper() == "Y":
            result = True
            break
        elif input_key.upper() == "N":
            result = False
            break
    return result


def input_labels(deuterated):
    text = ("-------\n"
            "Please specify amino acid labeling types with one letter codes\n"
            "Labeling_code X C A F N S T D  X F N D\n"
            "      Alpha-N 0 0 0 0 1 1 1 1  0 0 1 1\n"
            "           CO 0 1 0 1 0 1 0 1  0 1 0 1\n"
            "           CA 0 0 1 1 0 0 1 1  0 1 0 1\n"
            "   Deuterated 0 0 0 0 0 0 0 0  1 1 1 1\n"
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
            print("Please retry input")
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
            "5. COHNCA\n"
            "6. DQHNCA\n"
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
            print("Please retry input")
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


def input_name():
    repeat_input = True
    name = ''
    while repeat_input:
        input_line = input("------\nPlease enter the name of the NCS (Only letters or digits can be used): ")
        good_name = True
        for char in input_line:
            if not (char.isalpha() or char.isdigit()):
                print("Character '{}' can't be used in NCS name (not a digit or a letter)".format(char))
                good_name = False
                break
        if good_name:
            name = input_line
            print("Entered name is '{}'".format(name))
            break
        else:
            print("Please repeat the input")
            repeat_input = ask_to_continue_input()
    return name, ask_to_continue_input()


def print_description():
    description = ("This script is designed to create NMR Coding System (NCS) file\n"
                   "suitable for UUCSL package. To create a file you will need\n"
                   "to specify amino acid label types, are they deuterated or not,\n"
                   "NMR spectra and the name of the generated NCS\n"
                   )
    print(description)
    input("Press Enter to continue...")


def write_to_file(ncs):
    while True:
        input_key = input("Do you want to save this ncs to file? (Y/N) > ")
        if input_key.upper() == "Y":
            with open("{}.ncs".format(ncs.name), "w") as f:
                f.write(str(ncs))
                print("The '{}' NCS was successfully written to '{}.ncs'".format(ncs.name, ncs.name))
                break
        elif input_key.upper() == "N":
            break


def main():
    print_description()
    deuterated, continue_input = input_deuteration()
    if not continue_input:
        return
    labels, continue_input = input_labels(deuterated)
    if not continue_input:
        return
    spectra, continue_input = input_spectra()
    if not continue_input:
        return
    name, continue_input = input_name()
    if not continue_input:
        return
    ncs = NCS(name, spectra, labels)
    print(ncs)

    write_to_file(ncs)


if __name__ == "__main__":
    main()