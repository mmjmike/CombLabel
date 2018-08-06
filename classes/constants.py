class LabelType:

    def __init__(self, name, isotopes):
        self.name = name
        self.isotopes = isotopes


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

        else:
            print("Error")
            sys.stdout.flush()
            return 0


class NCS:

    NITRO_TYPES = ("N", "D", "S", "T")  # labeling types with 15N

    def __init__(self, name, spectra_list, label_types):
        self.name = name
        self.spec_list = spectra_list
        self.label_types = label_types
        self.label_dict = {}
        for label in self.label_types:
            self.label_dict.update({label.name: label})
        self.letters = ("a", "b", "c", "d", "e", "f")
        self.codes_dict = {}
        self.label_power = {}
        self.spectra_numbers = []
        self._make_coding_table()

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

    CARBON_TYPES = ("C", "D", "A", "S", "T", "F")  # labeling types with 13C

    HSQC = Spectrum("HSQC")
    HNCO = Spectrum("HNCO")
    HNCA = Spectrum("HNCA")
    HNCOCA = Spectrum("HNCOCA")
    DQHNCA = Spectrum("DQHNCA")
    COHNCA = Spectrum("COHNCA")
    HNCACO = Spectrum("HNCACO")

    basic_spectra = (HSQC, HNCO, HNCA, HNCOCA, COHNCA, DQHNCA, HNCACO)

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

    NC2       = NCS("NC2", [HSQC, HNCO],
                    [typeX, typeN, typeC])
    NC2nox    = NCS("NC2nox", [HSQC, HNCO],
                    [       typeN, typeC])
    NCD2      = NCS("NCD2", [HSQC, HNCO],
                    [typeX, typeN, typeC, typeD])
    NCD2nox   = NCS("NCD2nox", [HSQC, HNCO],
                    [       typeN, typeC, typeD])
    NCD4      = NCS("NCD4", [HSQC, HNCO, HNCA],
                    [typeX, typeN, typeC, typeD])
    NCD4nox   = NCS("NCD4nox", [HSQC, HNCO, HNCA],
                    [       typeN, typeC, typeD])
    NCD6      = NCS("NCD6", [HSQC, HNCO, HNCA, HNCOCA, DQHNCA],
                    [typeX, typeN, typeC, typeD])
    NCD6nox   = NCS("NCD6nox", [HSQC, HNCO, HNCA, HNCOCA, DQHNCA],
                    [       typeN, typeC, typeD])
    NCDF6      = NCS("NCDF6", [HSQC, HNCO, HNCA, HNCOCA, DQHNCA],
                    [typeX, typeN, typeC, typeD, typeF])
    NCDF6nox   = NCS("NCDF6nox", [HSQC, HNCO, HNCA, HNCOCA, DQHNCA],
                    [       typeN, typeC, typeD, typeF])
    NCDA8     = NCS("NCDA8", [HSQC, HNCO, HNCA, HNCOCA, COHNCA, DQHNCA],
                    [typeX, typeN, typeC, typeD,        typeA])
    NCDA8nox  = NCS("NCDA8nox", [HSQC, HNCO, HNCA, HNCOCA, COHNCA, DQHNCA],
                    [       typeN, typeC, typeD,        typeA])
    NCDAF8     = NCS("NCDAF8", [HSQC, HNCO, HNCA, HNCOCA, COHNCA, DQHNCA],
                    [typeX, typeN, typeC, typeD, typeF, typeA])
    NCDAF8nox  = NCS("NCDAF8nox", [HSQC, HNCO, HNCA, HNCOCA, COHNCA, DQHNCA],
                    [       typeN, typeC, typeD, typeF, typeA])
    TSF12     = NCS("TSF12", [HSQC, HNCO, HNCA, HNCOCA, COHNCA, DQHNCA, HNCACO],
                    [typeX, typeN, typeC, typeD, typeF, typeA, typeT, typeS])
    TSF12nox  = NCS("TSF12nox", [HSQC, HNCO, HNCA, HNCOCA, COHNCA, DQHNCA, HNCACO],
                    [       typeN, typeC, typeD, typeF, typeA, typeT, typeS])
    TF12      = NCS("TF12", [HSQC, HNCO, HNCA, HNCOCA, COHNCA, DQHNCA, HNCACO],
                    [typeX, typeN, typeC, typeD, typeF, typeA, typeT,      ])
    TSF12noxs = NCS("TSF12noxs", [HSQC, HNCO, HNCA, HNCOCA, COHNCA, DQHNCA, HNCACO],
                    [       typeN, typeC, typeD, typeF, typeA, typeT,      ])
    XND2      = NCS("XND2", [HSQC, HNCO],
                    [typeX, typeN, typeD])
    XND4      = NCS("XND4", [HSQC, HNCO, HNCA],
                    [typeX, typeN, typeD])

    LIST_OF_NCS = [NC2,    NC2nox, 
                   NCD2,   NCD2nox, 
                   NCD4,   NCD4nox,
                   NCD6,   NCD6nox,
                   NCDF6,  NCDF6nox,
                   NCDA8,  NCDA8nox,
                   NCDAF8, NCDAF8nox,
                   TSF12,  TSF12nox, TF12, TSF12noxs,
                   XND2, XND4]

    NCS_NAMES = [ncs.name for ncs in LIST_OF_NCS]

    def __init__(self):
        pass
