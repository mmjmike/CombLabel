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

        elif self.name == "COHNCA":
            return int(atom_list[3] and atom_list[1] and not atom_list[2])

        elif self.name == "DQHNCA":
            return int(atom_list[3] and atom_list[1] and atom_list[4])

        elif self.name == "HNCACO":
            return int(atom_list[3] and atom_list[4] and atom_list[5])


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
        self.letters = ("a", "b", "c", "d", "e", "f")
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

        return output


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
    CARBON_TYPES = ("C", "D", "A", "S", "T", "F")  # labeling types with 13C

    HSQC = Spectrum("HSQC")
    HNCO = Spectrum("HNCO")
    HNCA = Spectrum("HNCA")
    HNCOCA = Spectrum("HNCOCA")
    DQHNCA = Spectrum("DQHNCA")
    COHNCA = Spectrum("COHNCA")
    HNCACO = Spectrum("HNCACO")

    basic_spectra = (HSQC, HNCO, HNCA, HNCOCA, DQHNCA, COHNCA, HNCACO)

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
    TO_THREE_LETTER_CODE = {}
    for code3 in TO_ONE_LETTER_CODE:
        code1 = TO_ONE_LETTER_CODE[code3]
        TO_THREE_LETTER_CODE[code1] = code3

    RES_TYPES_THREE = (TO_ONE_LETTER_CODE[res] for res in RES_TYPES_LIST)

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
    NCDA8 = NCS("NCDA8", [HSQC, HNCO, HNCA, HNCOCA, COHNCA, DQHNCA],
                [typeX, typeN, typeC, typeD, typeA])
    TSF12 = NCS("TSF12", [HSQC, HNCO, HNCA, HNCOCA, COHNCA, DQHNCA, HNCACO],
                [typeX, typeN, typeC, typeD, typeA, typeT, typeS, typeF])
    TF12 = NCS("TF12", [HSQC, HNCO, HNCA, HNCOCA, COHNCA, DQHNCA, HNCACO],
                [typeX, typeN, typeC, typeD, typeA, typeT, typeF])
    XND2 = NCS("XND2", [HSQC, HNCO],
               [typeX, typeN, typeD])
    XND4 = NCS("XND4", [HSQC, HNCO, HNCA],
               [typeX, typeN, typeD])

    LIST_OF_NCS = [NC2, NCD2, NCD4, NCD6, NCDA8, TSF12, TF12, XND2, XND4]

    NCS_NAMES = [ncs.name for ncs in LIST_OF_NCS]

    def __init__(self):
        pass
