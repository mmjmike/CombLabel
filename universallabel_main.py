import copy

class Scheme:

    def __init__(self, name, ncs, samples, patterns):
        self.name = name
        self.ncs = ncs
        self.codes = set()
        self.patterns = patterns
        self.samples = samples
        self.good = self.check_codes()
        self.simplified = {}
        self.simplify()
        self.new_codes = set()

    def check_patterns(self, patterns):
        if not patterns:
            return False
        size = patterns[0].size
        for pattern in patterns:
            for label_type in pattern.string:
                if label_type not in self.ncs.label_types:
                    return False
            if pattern.size != size:
                return False
        return True

    def check_codes(self):
        self.codes = set()
        first = True
        for pattern_1 in self.patterns:
            for pattern_2 in self.patterns:
                code = self.ncs.calc_code(pattern_1, pattern_2)
                if first:
                    first = False
                    continue
                elif code in self.codes:
                    return False
                else:
                    self.codes.add(code)
        return True

    def add_pattern(self, new_pattern):
        if self.try_pattern(new_pattern):
            self.patterns.append(new_pattern)
            self.codes.update(self.new_codes)
            # if len(self.patterns) == 2:
            #     print(self.patterns)

    def try_pattern(self, new_pattern):
        if not self.good:
            return False
        self.new_codes = set()
        if new_pattern in self.patterns:
            return False
        for pattern in self.patterns:
            code_1 = self.ncs.calc_code(pattern, new_pattern)
            code_2 = self.ncs.calc_code(new_pattern, pattern)
            if code_1 in self.codes \
                    or code_2 in self.codes \
                    or code_1 == code_2 \
                    or code_1 in self.new_codes \
                    or code_2 in self.new_codes:
                return False
            else:
                self.new_codes.add(code_1)
                self.new_codes.add(code_2)
        self_code = self.ncs.calc_code(new_pattern, new_pattern)
        if self_code in self.new_codes or self_code in self.codes:
            return False
        else:
            self.new_codes.add(self_code)
        return True

    def add_pattern_list(self, pattern_list):
        for pattern in pattern_list:
            self.add_pattern(pattern)

    def simplify(self):
        for pattern in self.patterns:
            simple_pattern = self.simplify_pattern(pattern)
            if self.simplified != {} and simple_pattern in self.simplified:
                self.simplified[simple_pattern] += 1
            else:
                self.simplified.update({simple_pattern: 1})

    def simplify_pattern(self, pattern):
        simple_form = [0 for _ in range(len(self.ncs.TYPES))]
        for label in pattern:
            for i in range(len(self.ncs.TYPES)):
                if self.ncs.TYPES[i] == label:
                    simple_form[i] += 1
                    continue
        result = "".join([str(a) for a in simple_form])
        return result

    def equal_to(self, scheme):
        return self.simplified == scheme.simplified

    def copy(self):
        return Scheme(copy.copy(self.name), copy.copy(self.ncs),
                      copy.copy(self.samples), copy.copy(self.patterns))

    def direct_product(self, scheme):
        new_patterns = []
        for pattern_1 in self.patterns:
            for pattern_2 in scheme.patterns:
                new_patterns.append(pattern_1 + pattern_2)
        new_name = self.name + "_X_" + scheme.name
        new_scheme = Scheme(new_name, self.ncs, new_patterns)
        return new_scheme

    def clear(self):
        self.codes = set()
        self.patterns = []
        self.size = 0
        self.good = True
        self.simplified = {}
        self.simplify()
        self.new_codes = set()

    def output(self):
        output = "\n"
        for pattern in self.patterns:
            output += "\n" + pattern
        print(output)

    def size(self):
        return len(self.patterns)

    def __str__(self):
        output = ""
        for pattern in self.patterns:
            output += "\n" + pattern
        return output


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
            return 0


class NCS:

    TYPES = ("X", "N", "C", "D", "A", "S", "T", "F")
    NITRO_TYPES = ("N", "D", "S", "T")  # labeling types with 15N
    CARBON_TYPES = ("C", "D", "A", "S", "T", "F")  # labeling types with 13C

    def __init__(self, spectra_list, label_types):
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


class LabelType:

    def __init__(self, name, isotopes):
        self.name = name
        self.isotopes = isotopes


class Task:

    def __init__(self, types, samples, ncs, min_depth, target_depth):
        self.types = types
        self.samples = samples
        self.ncs = ncs
        self.min_depth = min_depth
        self.target_depth = target_depth
        self.scheme = Scheme("1", self.ncs, samples, [])
        self.patterns = [list(self.generate_patterns(self.samples))]
        self.depth = 0
        self.counter = [0]
        self.back_up_schemes = []
        self.result = {}
        # print(self.patterns)

    def zg(self):
        while True:
            patterns = self.patterns[self.depth]
            if self.depth == 0 and self.counter[0] + 1 == len(patterns):
                break
            self.back_up_schemes.append(self.scheme.copy())
            self.scheme.add_pattern(patterns[self.counter[self.depth]])
            next_patterns = []
            start_point = 1 + self.counter[self.depth]
            patterns_left = len(patterns) - start_point
            # print(patterns_left, self.depth)
            if patterns_left == 0:
                if len(self.scheme.patterns) >= self.min_depth:
                    self.save_result()
                self.depth -= 1
                self.patterns.pop()
                self.counter.pop()
                self.counter[-1] += 1
                self.back_up_schemes.pop()
                self.scheme = self.back_up_schemes[-1].copy()
                self.back_up_schemes.pop()
                continue
            for i in range(patterns_left):
                if self.scheme.try_pattern(patterns[i+start_point]):
                    next_patterns.append(patterns[i+start_point])
            if next_patterns == []:

                self.scheme = self.back_up_schemes[self.depth].copy()
                self.back_up_schemes.pop()
                self.counter[self.depth] += 1
            else:
                self.patterns.append(next_patterns)
                # print(next_patterns, self.depth)
                self.counter.append(0)
                self.depth += 1
        self.write_result()

    def generate_patterns(self, samples, top=True):
        if samples == 0:
            new_set = set()
            new_set.add("")
            return new_set
        current_set = self.generate_patterns(samples - 1, top=False)
        new_set = set()
        for item in current_set:
            for option in self.types:
                new_pattern = item + option.name
                if top:
                    if self.ncs.check_power(new_pattern, self.min_depth):
                        new_set.add(new_pattern)
                else:
                    new_set.add(new_pattern)
        return new_set

    def save_result(self):
        depth_of_scheme = self.scheme.size()
        if depth_of_scheme in self.result:
            for scheme in self.result[depth_of_scheme]:
                if not self.scheme.equal_to(scheme):
                    self.result[depth_of_scheme].append(self.scheme.copy())
        else:
            self.result.update({depth_of_scheme: [self.scheme.copy()]})

    def write_result(self):
        for depth, schemes in self.result:
            for scheme in schemes:
                print(depth, "\n", scheme)


typeX = LabelType("X", "000")
typeN = LabelType("N", "100")
typeC = LabelType("C", "001")
typeD = LabelType("D", "111")
typeA = LabelType("A", "010")
typeT = LabelType("T", "110")
typeS = LabelType("S", "101")
typeF = LabelType("F", "011")
basic_types = (typeX, typeN, typeC, typeD, typeA, typeT, typeS, typeF)

HSQC = Spectrum("HSQC")
HNCO = Spectrum("HNCO")
HNCA = Spectrum("HNCA")
HNCOCA = Spectrum("HNCOCA")
DQHNCA = Spectrum("DQHNCA")
COHNCA = Spectrum("COHNCA")
HNCACO = Spectrum("HNCACO")

basic_spectra = (HSQC, HNCO, HNCA, HNCOCA, COHNCA, DQHNCA, HNCACO)

NC2 = NCS([HSQC, HNCO], [typeX, typeN, typeC])

task = Task([typeX, typeN, typeC], 3, NC2, 3, 4)

def main():
    task.zg()


if __name__ == "__main__":
    main()