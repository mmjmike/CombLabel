import re


class LabelType:

    def __init__(self, name, isotopes):
        self.name = name
        self.isotopes = isotopes

    def __str__(self):
        return self.name


class Spectrum:

    def __init__(self, name):
        self.name = name

    def has_signal(self, label_type_1, label_type_2):
        atom_list = [int(i) for i in (label_type_1.isotopes + label_type_2.isotopes)]
        if self.name == "HSQC":
            return int(atom_list[3])

        elif self.name == "HNCO":
            return int(atom_list[3] and atom_list[2])

        elif self.name == "HNCA":
            return int(atom_list[3] and (atom_list[1] or atom_list[4]))

        elif self.name == "HNCOCA":
            return int(atom_list[3] and atom_list[2] and atom_list[1])

        elif self.name == "COfHNCA":
            return int(atom_list[3] and atom_list[1] and not atom_list[2])

        elif self.name == "DQHNCA":
            return int(atom_list[3] and atom_list[1] and atom_list[4])

        elif self.name == "HNCACO":
            return int(atom_list[3] and atom_list[4] and atom_list[5])

        elif self.name == "HNCAfCO":
        # Hypothetic spectrum, to be devepoled
        # HN(CA-filtered)CO
            return int(atom_list[3] and not atom_list[4] and atom_list[5])


class NCS:

    NITRO_TYPES = ("N", "D", "S", "T")  # labeling types with 15N

    def __init__(self, name, spectra_list, label_types, deuterated=False):
        self.name = name
        self.deuterated = deuterated
        self.spec_list = spectra_list
        self.label_types = label_types
        self.label_dict = {}
        for label in self.label_types:
            self.label_dict.update({label.name: label})
        self.letters = ("a", "b", "c", "d", "e", "f", "g", "h", "i", "j")
        self.codes_dict = {}
        self.label_power = {}
        self.spectra_numbers = []
        self.vectors = [[0 for _ in self.spec_list]]
        self._make_coding_table()

    @property
    def max_code_value(self):
        return len(self.vectors) - 1

    def _make_coding_table(self):

        codes_table = [[0 for i in range(len(self.label_types))] for j in range(len(self.label_types))]
        self.vectors = [[0 for _ in self.spec_list]]
        for i in range(len(self.label_types)):
            label_2 = self.label_types[i]
            for j in range(len(self.label_types)):
                label_1 = self.label_types[j]
                vector = []
                for spectrum in self.spec_list:
                    vector.append(spectrum.has_signal(label_1, label_2))
                if vector in self.vectors:
                    code = self.vectors.index(vector)
                else:
                    code = len(self.vectors)
                    self.vectors.append(vector)
                if code > 9:
                    codes_table[j][i] = self.letters[code-10]
                else:
                    codes_table[j][i] = str(code)
        for i in range(len(self.label_types)):
            result = []
            for row in codes_table:
                result.append(row[i])
            power = len(set(result))
            self.label_power[self.label_types[i]] = power
            if power == 1 and self.label_types[i] in self.NITRO_TYPES:
                raise ValueError
        for i in range(len(self.label_types)):
            label_1 = self.label_types[i]
            subdict = {}
            for j in range(len(self.label_types)):
                label_2 = self.label_types[j]
                subdict[label_2] = codes_table[i][j]
            self.codes_dict[label_1] = subdict

    def calc_code(self, pattern_1, pattern_2):
        return ''.join([self.codes_dict[self.label_dict[pattern_1[i]]][self.label_dict[pattern_2[i]]] for i in range(len(pattern_1))])

    def check_power(self, pattern, min_power):
        power = 1
        for label in pattern:
            power *= self.label_power[self.label_dict[label]]
        return power >= min_power

    def __str__(self):
        output = "[NCS = {}]\n".format(self.name)
        output += "[Deuterated = {}]\n".format(self.deuterated)
        output += "[Labels = {}]\n".format(",".join([label.name for label in self.label_types]))
        output += "[Spectra = {}]\n".format(",".join([spectrum.name for spectrum in self.spec_list]))
        output += "\n\n" + "#" * 50 + "\n"
        output += "# Spectrum code for each labeling pair \n#\n"
        output += "# One-letter codes in the headers of columns"
        output += "# and rows are the labeling types \n"
        output += "# Don't confuse with one-letter amino acid codes\n\n"
        output += "[code_pairs]\n "
        for type in self.label_types:
            output += ", " + type.name
        output += "\n"
        for type_1 in self.label_types:
            output += type_1.name
            for type_2 in self.label_types:
                output += ", " + self.codes_dict[type_1][type_2]
            output += "\n"
        output += "\n"

        output += "\n\n"+"#"*50+"\n"
        output += "# Spectrum codes table\n#\n"
        output += "# The spectrum code is in the first column\n"
        output += "# Flag of the peak presence (0 or 1) for each spectrum\n\n"
        output += "[codes]\n"
        output += "Code"
        for spectrum in self.spec_list:
            output += "," + "{:>7}".format(spectrum.name)
        output += "\n"
        for i in range(len(self.vectors)):
            if i < 10:
                code = "{:>4}".format(i)
            else:
                code = "{:>4}".format(self.letters[i-10])
            output += code + ","
            output += ",".join(["{:>7}".format(item) for item in self.vectors[i]])
            if i+1 < len(self.vectors):
                output += "\n"
        output += "\n\n"

        return output

class ELB:

    def __init__(self, patterns, ncs_name, deuterated=False):
        self.patterns = patterns
        self.ncs_name = ncs_name
        self.deuterated = deuterated
        self.simplified = {}
        self.simplified_str = ""
        self.simplify()

    def __str__(self):
        return "\n".join(self.patterns)

    @property
    def samples(self):
        if self.patterns:
            return len(self.patterns[0])
        else:
            return 0

    def full_str(self):
        return "[ELB samples = {} patterns = {}]\n".format(self.samples, len(self.patterns)) \
                 + str(self)

    def simplify(self):
        self.simplified, self.simplified_str = \
            Pattern.simplify_list_of_patterns(self.patterns)

    def __mul__(self, other):
        new_patterns = []
        for pattern_1 in self.patterns:
            for pattern_2 in other.patterns:
                new_patterns.append(pattern_1 + pattern_2)
        return ELB(new_patterns, self.ncs_name, self.deuterated)

    def __eq__(self, scheme):
        return self.simplified == scheme.simplified

    def __hash__(self):
        return(hash(self.simplified_str))

    def sort(self):
        for i in range(len(self.patterns)-1):
            for j in range(len(self.patterns)-1-i):
                if Pattern.pattern_bigger(self.patterns[i], self.patterns[i+j+1]):
                    temp_pattern = self.patterns[i]
                    self.patterns[i] = self.patterns[i+j+1]
                    self.patterns[i+j+1] = temp_pattern

    def is_subset_of(self, other_simple):
        for pattern in self.simplified:
            if pattern not in other_simple:
                return False
            if self.simplified[pattern] > other_simple[pattern]:
                return False
        return True



class ELB_to_be_removed:

    def __init__(self, patterns, ncs_name, deuterated=False):
        self.patterns = patterns
        self.ncs_name = ncs_name
        self.deuterated = deuterated
        self.simplified = {}
        self.simplify()

    def __str__(self):
        return "\n".join(self.patterns)

    @property
    def samples(self):
        if self.patterns:
            return len(self.patterns[0])
        else:
            return 0

    def full_str(self):
        return "[ELB samples = {} patterns = {}]\n".format(self.samples, len(self.patterns)) \
                 + str(self)

    def __mul__(self, other):
        new_patterns = []
        for pattern_1 in self.patterns:
            for pattern_2 in other.patterns:
                new_patterns.append(pattern_1 + pattern_2)
        return ELB(new_patterns, self.ncs_name, self.deuterated)

    def __eq__(self, scheme):
        return self.simplified == scheme.simplified

    def simplify(self):
        self.simplified = {}
        for pattern in self.patterns:
            simple_pattern = Pattern.simplify_pattern(pattern)
            if self.simplified != {} and simple_pattern in self.simplified:
                self.simplified[simple_pattern] += 1
            else:
                self.simplified.update({simple_pattern: 1})

    def sort(self):
        for i in range(len(self.patterns)-1):
            for j in range(len(self.patterns)-1-i):
                if Pattern.pattern_bigger(self.patterns[i], self.patterns[i+j+1]):
                    temp_pattern = self.patterns[i]
                    self.patterns[i] = self.patterns[i+j+1]
                    self.patterns[i+j+1] = temp_pattern

    def is_subset_of(self, other_simple):
        for pattern in self.simplified:
            if pattern not in other_simple:
                return False
            if self.simplified[pattern] > other_simple[pattern]:
                return False
        return True


class Constants:

    typeX = LabelType("X", "000")
    typeN = LabelType("N", "100")
    typeC = LabelType("C", "001")
    typeD = LabelType("D", "111")
    typeA = LabelType("A", "010")
    typeT = LabelType("T", "110")
    typeS = LabelType("S", "101")
    typeF = LabelType("F", "011")
    BASIC_TYPES = (typeX, typeN, typeC, typeD, typeA, typeT, typeS, typeF)
    TYPES_NAMES = [label_type.name for label_type in BASIC_TYPES]

    TYPES = ("X", "N", "C", "D", "A", "T", "S", "F")
    DEUTERATED_TYPES = ("X", "N", "F", "D")
    CARBON_TYPES = [label_type for label_type in BASIC_TYPES if label_type.isotopes[1] == "1" or label_type.isotopes[2] == "1"]
    NITRO_TYPES = [label_type for label_type in BASIC_TYPES if label_type.isotopes[0] == "1"]

    HSQC =    Spectrum("HSQC")
    HNCO =    Spectrum("HNCO")
    HNCA =    Spectrum("HNCA")
    HNCOCA =  Spectrum("HNCOCA")
    DQHNCA =  Spectrum("DQHNCA")
    COfHNCA = Spectrum("COfHNCA")
    HNCACO =  Spectrum("HNCACO")
    HNCAfCO = Spectrum("HNCAfCO")

    basic_spectra = (HSQC, HNCO, HNCA, HNCOCA, DQHNCA, COfHNCA, HNCACO, HNCAfCO)
    SPECTRA_NAMES = [spectrum.name for spectrum in basic_spectra]
    SPECTRA_NAMES_UPPER = [spectrum.name.upper() for spectrum in basic_spectra]
  

    RES_TYPES_LIST = ("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                      "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

    TO_ONE_LETTER_CODE = {
        "Ala": "A", "Cys": "C",
        "Asp": "D", "Glu": "E",
        "Phe": "F", "Gly": "G",
        "His": "H", "Ile": "I",
        "Lys": "K", "Leu": "L",
        "Met": "M", "Asn": "N",
        "Pro": "P", "Gln": "Q",
        "Arg": "R", "Ser": "S",
        "Thr": "T", "Val": "V",
        "Trp": "W", "Tyr": "Y"
    }



    RES_TYPES_THREE = []
    TO_THREE_LETTER_CODE = {}
    for code3 in TO_ONE_LETTER_CODE:
          code1 = TO_ONE_LETTER_CODE[code3]
          TO_THREE_LETTER_CODE[code1] = code3
          RES_TYPES_THREE.append(code3)


    PROLINE_SUBSTITUTE = {
        "X": "X",
        "N": "X",
        "C": "C",
        "D": "F",
        "F": "F",
        "T": "A",
        "A": "A",
        "S": "C"
    }

    NC2 = NCS("NC2", [HSQC, HNCO],
              [typeX, typeN, typeC])
    NCD2 = NCS("NCD2", [HSQC, HNCO],
               [typeX, typeN, typeC, typeD])
    NCD4 = NCS("NCD4", [HSQC, HNCO, HNCA],
               [typeX, typeN, typeC, typeD])
    NCD6 = NCS("NCD6", [HSQC, HNCO, HNCA, HNCOCA, DQHNCA],
               [typeX, typeN, typeC, typeD])
    NCDA8 = NCS("NCDA8", [HSQC, HNCO, HNCA, HNCOCA, COfHNCA, DQHNCA],
                [typeX, typeN, typeC, typeD, typeA])
    TSF12 = NCS("TSF12", [HSQC, HNCO, HNCA, HNCOCA, COfHNCA, DQHNCA, HNCACO],
                [typeX, typeN, typeC, typeD, typeA, typeT, typeS, typeF])
    TF12 = NCS("TF12", [HSQC, HNCO, HNCA, HNCOCA, COfHNCA, DQHNCA, HNCACO],
                [typeX, typeN, typeC, typeD, typeA, typeT, typeF])
    XND2 = NCS("XND2", [HSQC, HNCO],
               [typeX, typeN, typeD])
    XND4 = NCS("XND4", [HSQC, HNCO, HNCA],
               [typeX, typeN, typeD])

    LIST_OF_NCS = [NC2, NCD2, NCD4, NCD6, NCDA8, TSF12, TF12, XND2, XND4]

    NCS_NAMES = [ncs.name for ncs in LIST_OF_NCS]

    ncs_re = re.compile('\\[\\s*NCS\\s*=\\s*(\\w+)\\s*\\]')
    deuterated_re = re.compile('\\[\\s*Deuterated\\s*=\\s*([A-Za-z]+)\\s*\\]')
    elb_re = re.compile('\\[\\s*ELB\\s+samples\\s*=\\s*(\\d+)\\s+patterns\\s*=\\s*(\\d+)\\s*\\]')

    labels_re = re.compile('\\[\\s*Labels\\s*=\\s*([A-Za-z ,]+)\\s*\\]')
    spectra_re = re.compile('\\[\\s*Spectra\\s*=\\s*([A-Za-z0-9 ,]+)\\s*\\]')
    solution_re = re.compile('\\[\\s*solution\\s*\\]')


class PatternClass:

    def pattern_bigger(self, pattern1, pattern2):
        for i in range(len(pattern1)):
            type1 = pattern1[i]
            type2 = pattern2[i]
            if Constants.TYPES.index(type1) > Constants.TYPES.index(type2):
                return True
            if Constants.TYPES.index(type1) < Constants.TYPES.index(type2):
                return False
        return True

    def simplify_pattern(self, pattern):
        simple_form = [0 for _ in range(len(Constants.TYPES))]
        for label in pattern:
           for i in range(len(Constants.TYPES)):
               if Constants.TYPES[i] == label:
                   simple_form[i] += 1
                   continue
        letters = list(map(chr, range(97, 123)))
        result = "".join([str(a) if a < 10 else letters[a-10] for a in simple_form])
        return result

    def simplify_list_of_patterns(self, list_of_patterns):
        simplified = {}
        for pattern in list_of_patterns:
            simple_pattern = self.simplify_pattern(pattern)
            if simplified != {} and simple_pattern in simplified:
                simplified[simple_pattern] += 1
            else:
                simplified.update({simple_pattern: 1})
        out = []
        for simple_pattern in sorted(simplified.keys()):
            out.append(simple_pattern + ":" + str(simplified[simple_pattern]))
        simplified_str = ",".join(out)
        return simplified, simplified_str

    def first_scheme_subset(self, scheme_1, scheme_2):
        for pattern in scheme_1:
            if pattern not in scheme_2:
                return False
            if scheme_1[pattern] > scheme_2[pattern]:
                return False
        return True

    def count_type_in_list_of_patterns(self, patterns, label_type):
        simplified, simplified_str = self.simplify_list_of_patterns(patterns)
        index_of_t = self.index_of_type(label_type)
        return self.count_type_in_list_of_simplified(simplified, index_of_t)

    def count_type_in_list_of_simplified(self, simplified, index_of_type):
        count_type = 0
        count_all  = 0
        for simple_pattern, pattern_count in simplified.items():
            has_t = 0
            if int(simple_pattern[index_of_type]):
                has_t = 1
            count_type += has_t * pattern_count
            count_all  += pattern_count
        return count_type, count_all - count_type

    def index_of_type(self, label_type):
        index_of_t = -1
        try:
            index_of_t = Constants.BASIC_TYPES.index(label_type)
        except:
            print("Internal error: labeling type \"{}\" is not found in Constants.BASIC_TYPES. Exiting.".format(
                label_type.name))
            exit(-1)
        return index_of_t

    def __init__(self):
        pass


def auto_name_ncs(ncs):
    auto_name = ''
    if ncs.deuterated:
        auto_name += "2H-"
    flag_X = False
    for label_type in ncs.label_types:
        if label_type.name == "X":
            flag_X = True
        if not ncs.deuterated and label_type.name == "X":
            continue
        auto_name += label_type.name
    auto_name += str(ncs.max_code_value)
    if not ncs.deuterated and not flag_X:
        auto_name += "noX"
    return auto_name

Pattern = PatternClass()
