#!/usr/bin/python3

import copy
from scipy.optimize import linprog
from cl_errors import errors as err


PRICES_FILE = "BEST_PRICE_1H.csv"
FULL_OUTPUT = True
JOB_NAME = "TSF12-2"
OUTPUT_FILE = JOB_NAME + ".txt"
WRITE_TO_FILE = True
WRITE_TO_CONSOLE = True


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
        self.precursor_schemes = []
        self.precursor_block_types = []

    # def check_patterns(self, patterns):
    #     if not patterns:
    #         return False
    #     size = patterns[0].size
    #     for pattern in patterns:
    #         for label_type in pattern.string:
    #             if label_type not in self.ncs.label_types:
    #                 return False
    #         if pattern.size != size:
    #             return False
    #     return True

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

    def sort(self):
        for i in range(len(self.patterns)-1):
            for j in range(len(self.patterns)-1-i):
                if self.pattern_bigger(self.patterns[i], self.patterns[i+j+1]):
                    temp_pattern = self.patterns[i]
                    self.patterns[i] = self.patterns[i+j+1]
                    self.patterns[i+j+1] = temp_pattern

    def pattern_bigger(self, pattern1, pattern2):
        for i in range(len(pattern1)):
            type1 = pattern1[i]
            type2 = pattern2[i]
            if self.ncs.TYPES.index(type1) > self.ncs.TYPES.index(type2):
                return True
            if self.ncs.TYPES.index(type1) < self.ncs.TYPES.index(type2):
                return False
        return True

    def add_pattern(self, new_pattern):
        if self.try_pattern(new_pattern):
            self.patterns.append(new_pattern)
            self.codes.update(self.new_codes)
            self.simplify()
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
        self.simplified = {}
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
        # print(scheme.simplified)
        # print(self.simplified)
        a = (self.simplified == scheme.simplified)
        # print(a)
        return a

    def copy(self):
        return Scheme(copy.copy(self.name), copy.copy(self.ncs),
                      copy.copy(self.samples), copy.copy(self.patterns))

    def direct_product(self, scheme):
        new_patterns = []
        for pattern_1 in self.patterns:
            for pattern_2 in scheme.patterns:
                new_patterns.append(pattern_1 + pattern_2)
        new_name = self.name + "_X_" + scheme.name
        samples = len(new_patterns[0])
        new_scheme = Scheme(new_name, self.ncs, samples, new_patterns)
        return new_scheme

    def clear(self):
        self.codes = set()
        self.patterns = []
        # self.size = 0
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


class BlockFinder:

    def __init__(self, types, samples, ncs, min_depth):
        self.types = types
        self.samples = samples
        self.ncs = ncs
        self.min_depth = min_depth
        self.scheme = Scheme("1", self.ncs, samples, [])
        self.patterns = [list(self.generate_patterns(self.samples))]
        self.depth = 0
        self.counter = [0]
        self.back_up_schemes = []
        self.result = {}
        self.results_found = 0
        self.output = ''
        # print(self.patterns)

    def find(self):
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
                    # self.scheme.output()
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
        self.scheme.sort()

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
            equal = False
            for scheme in self.result[depth_of_scheme]:
                if self.scheme.equal_to(scheme):
                    equal = True
            if not equal:
                self.result[depth_of_scheme].append(self.scheme.copy())
                self.results_found += 1
                # print ("Elementary blocks found: ", self.results_found)

                    # print("Saved:")
                    # print(self.scheme.simplified)
                    # print(scheme.simplified)
                    #
                    # print(self.scheme)
        else:
            self.result.update({depth_of_scheme: [self.scheme.copy()]})

    def write_result(self):
        self.output = ''
        for depth in self.result:
            for scheme in self.result[depth]:
                # print("Output:")
                scheme.sort()
                self.output += '{} {}\n'.format(depth, scheme)


class PriceOptimizer:

    def __init__(self, ncs, price_filename):
        self.prices = {}
        self.ncs = ncs
        self.price_file = price_filename
        self.read_prices()

    def minimize_price(self, scheme, aa_list):
        patterns = []
        b_eq = [1 for _ in aa_list]
        b_ub = []
        A_eq = [[] for _ in range(len(aa_list))]
        price_sheet = []
        c = []
        x_coords = []
        for key in scheme.simplified:
            patterns.append(key)
            b_ub.append(scheme.simplified[key])
        A_ub = [[] for i in range(len(patterns))]
        for i in range(len(aa_list)):
            aa = aa_list[i]
            aa_prices = []
            for j in range(len(patterns)):
                pattern = patterns[j]
                price = self.calc_price(aa, pattern)
                aa_prices.append(price)
                if price > 0:
                    for k in range(len(aa_list)):
                        A_eq[k].append(0)
                    for k in range(len(patterns)):
                        A_ub[k].append(0)
                    A_eq[i][-1] = 1
                    A_ub[j][-1] = 1
                    c.append(price)
                    x_coords.append((i, j))
            price_sheet.append(aa_prices)

        result = linprog(c, A_ub, b_ub, A_eq, b_eq)
        # print(patterns)
        return [result, patterns, price_sheet]

    def calc_price(self, aa, pattern):
        price = 0
        for i in range(len(pattern)):
            label = pattern[i]
            number = int(label)
            type = self.ncs.TYPES[i]
            if aa == "P":
                if type == "N":
                    type = "X"
                if type == "D":
                    type = "F"
            curr_price = self.prices[aa][type] * number
            if curr_price < 0:
                return -1
            price += curr_price
        return price

    def read_prices(self):
        lines = self.get_lines_from_file(self.price_file)
        d = []
        self.prices = {}
        for line in lines:
            if line[0] != "#" and line != "":
                s = [x.strip() for x in line.split(",")]
                d.append(s)
        self.price_label_types = d[0][1:]
        # for label in self.ncs.TYPES:
        #     if label.name not in self.price_label_types:
        #         raise err.ReadConfigError(msg="ERROR! Prices are not specified for '" + label.name + "' label type")

        # ADD check if file format is correct
        # try:
        for i in range(len(d) - 1):
            curr_dict = {}
            for j in range(len(self.price_label_types)):
                try:
                    price = float(d[i + 1][j + 1])
                except ValueError:
                    raise err.ReadConfigError(
                        msg="ERROR! Price must be set in digits\nPlease check price file (row {}; col {})".format(
                            i + 2, j + 2))
                curr_dict.update({self.price_label_types[j]: price})
            residue_type = d[i + 1][0]
            self.prices.update({residue_type: curr_dict})
        # except IndexError:
        #     raise err.ReadConfigError(
        #         msg="ERROR in price file.\nPlease check the length of each row.\nUse comma as separator")

    def get_lines_from_file(self, filename):
        try:
            f = open(filename, 'r', encoding='UTF-8')
            lines = f.readlines()
            f.close()
        except IOError:
            raise err.ReadLinesError(filename)
        new_lines = []
        for line in lines:
            curr_line = line.rstrip()
            if curr_line:
                new_lines.append(curr_line)
        return new_lines


class Task:

    def __init__(self, ncs, aa_list, max_samples, max_block_size):
        self.ncs = ncs
        self.residues = aa_list
        self.aa_number = len(aa_list)
        self.max_samples = max_samples
        self.max_block_size = max_block_size
        self.results = {}
        self.types = []
        self.products = []
        self.all_schemes = []
        self.all_blocks = []
        self.product_schemes = 0
        self.first_output = True

    def find_scheme(self):
        self.find_blocks()
        self.find_products()
        checked_schemes = []
        checked_number = 0
        scheme_found = False
        price_optimizer = PriceOptimizer(self.ncs, PRICES_FILE)
        for product in self.products:
            self.all_schemes = []
            self.obtain_scheme(product, [Scheme("", self.ncs, 0, [""])], [[]])
            self.output("\nChecking product: {}".format(product))
            self.output("total schemes for this product: {}".format(len(self.all_schemes)))
            for i in range(len(self.all_schemes)):
                scheme = self.all_schemes[i]
                scheme.sort()
                blocks = self.all_blocks[i]
                self.output("\n{} scheme. Consists of following blocks: ".format(checked_number+1))
                for block in blocks:
                    self.output(str(block))
                self.output("\nScheme itself:")
                self.output(str(scheme))
                if scheme.simplified not in checked_schemes:
                    optimal_scheme = price_optimizer.minimize_price(scheme, self.residues)
                    # if impossible to use due to lack in the stock, skip
                    if not optimal_scheme[0].success:
                        self.output("Scheme can not be optimized")

                    elif checked_schemes == []:
                        best_scheme = optimal_scheme
                        best_scheme_patterns = copy.copy(scheme)
                        best_blocks = blocks
                        scheme_found = True
                    elif optimal_scheme[0].fun < best_scheme[0].fun:
                        best_scheme = optimal_scheme
                        best_scheme_patterns = copy.copy(scheme)
                        best_blocks = blocks
                    checked_schemes.append(scheme.simplified)
                else:
                    self.output("Equivalent scheme was already checked")
                checked_number += 1

                self.output("{}/{} schemes checked".format(checked_number,
                                                           self.product_schemes))
                if scheme_found:
                    self.output("Best price: {}".format(best_scheme[0].fun))

        if scheme_found:
            best_scheme_patterns.sort()
            label_dict = self.get_scheme_from_linprog(best_scheme, best_scheme_patterns)
            samples = best_scheme_patterns.samples
            self.output_best_scheme(label_dict, samples, best_blocks)
        else:
            self.output("No schemes were found...")

    def get_scheme_from_linprog(self, best_scheme, scheme_patterns):
        simple_patterns = best_scheme[1]
        simple_patterns_num = len(simple_patterns)
        result = best_scheme[0]
        price_sheet = best_scheme[2]
        x = result.x
        scheme = {}
        patterns = list(scheme_patterns.patterns)
        skip = 0
        for i in range(self.aa_number):
            res = self.residues[i]
            for j in range(simple_patterns_num):
                if price_sheet[i][j] < 0:
                    skip += 1
                    continue
                if x[i*simple_patterns_num+j-skip]:
                    simple_pattern = simple_patterns[j]
                    for k in range(len(patterns)):
                        pattern = patterns[k]
                        if scheme_patterns.simplify_pattern(pattern) == simple_pattern:
                            scheme[res] = ",".join(list(pattern))
                            patterns.pop(k)
                            break
        return scheme

    def output_best_scheme(self, scheme, samples, blocks):
        output = '\nBest scheme:\n[solution]\nRes'
        for i in range(samples):
            output += ",S{}".format(i + 1)
        output += "\n"
        for res in self.residues:
            output += "{}: {}\n".format(res, scheme[res])
        output += "\nBlocks used for this scheme"
        for block in blocks:
            output += str(block) + "\n"
        self.output(output)

    def obtain_scheme(self, product, schemes, block_list, depth=0):
        if depth >= len(product):
            self.all_schemes = schemes
            self.all_blocks = block_list
            return
        new_schemes = []
        new_block_list = []
        for scheme_1 in self.results[product[depth][1]]:
            for scheme_2 in schemes:
                new_schemes.append(scheme_2.direct_product(scheme_1))
        for scheme_1 in self.results[product[depth][1]]:
            for block in block_list:
                new_block_list.append(block + [scheme_1])
        self.obtain_scheme(product, new_schemes, new_block_list, depth=(depth+1))

    def add_block(self, curr_list_of_types, curr_size, curr_samples):
        for type in self.types:
            if curr_list_of_types and not self.types_in_order(curr_list_of_types[-1], type):
                continue
            size = curr_size * type[1]
            samples = curr_samples + type[0]
            if samples > self.max_samples:
                continue
            new_list_of_types = copy.copy(curr_list_of_types)
            new_list_of_types.append(type)
            if size >= self.aa_number:
                self.products.append(new_list_of_types)
            else:
                self.add_block(new_list_of_types, size, samples)

    def types_in_order(self, type_1, type_2):
        if type_2[0] > type_1[0]:
            return True
        elif type_2[0] == type_1[0]:
            if type_2[1] >= type_1[1]:
                return True
        return False

    def find_products(self):
        self.output("Calculating all possible product types")
        self.products = []
        self.add_block([], 1, 0)

        self.product_schemes = 0
        for list_of_types in self.products:
            number_of_schemes = 1
            for type in list_of_types:
                new_number = len(self.results[type[1]])
                number_of_schemes *= new_number
            self.product_schemes += number_of_schemes
        self.output("{} product types calculated".format(len(self.products)))
        self.output("{} total labeling schemes to check".format(self.product_schemes))

    def find_blocks(self):
        self.output("Search for elementary blocks with max {} samples".format(self.max_block_size))
        min_depth = 1
        self.results = {}
        self.types = []
        blocks_total = 0
        for i in range(self.max_block_size):
            self.output("Search in {} samples:".format(i+1))
            block_finder = BlockFinder(self.ncs.label_types, i+1, self.ncs, min_depth)
            block_finder.find()
            result = block_finder.result
            output = block_finder.output
            self.output(output)
            self.results.update(result)
            for key in result:
                self.types.append((i+1, key))
                blocks_total += len(result[key])
            # print(result)
            if result:
                min_depth = max(result) + 1
            else:
                min_depth += 1
        self.output("{} elementary blocks found\n".format(blocks_total))

    def print_blocks(self):
        output = ''
        for key in self.results:
            output += str(key) + "\n"
            output += self.results[key] + "\n"
        self.output(output)

    def output(self, text):
        mode = 'a'
        if self.first_output:
            mode = 'w'
            self.first_output = False
        if WRITE_TO_FILE:
            self.write_to_file(OUTPUT_FILE, text + "\n", mode=mode)
        if WRITE_TO_CONSOLE:
            print(text)

    def write_to_file(self, filename, output, mode='w'):
        f = open(filename, mode)
        f.write(output)
        f.flush()
        f.close()


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

NC2 = NCS([HSQC, HNCO],
          [typeX, typeN, typeC])
NCD2 = NCS([HSQC, HNCO],
           [typeX, typeN, typeC, typeD])
NCD4 = NCS([HSQC, HNCO, HNCA],
           [typeX, typeN, typeC, typeD])
NCD6 = NCS([HSQC, HNCO, HNCA, HNCOCA, DQHNCA],
           [typeX, typeN, typeC, typeD])
NCDA8 = NCS([HSQC, HNCO, HNCA, HNCOCA, COHNCA, DQHNCA],
            [typeX, typeN, typeC, typeD, typeA])
TSF12 = NCS([HSQC, HNCO, HNCA, HNCOCA, COHNCA, DQHNCA, HNCACO],
            [typeX, typeN, typeC, typeD, typeA, typeT, typeS, typeF])

RES_TYPES_LIST = ("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                  "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

task1 = Task(NC2, RES_TYPES_LIST, 9, 5)
task2 = Task(NCD2, RES_TYPES_LIST, 7, 4)
task3 = Task(NCD4, RES_TYPES_LIST, 6, 3)
task4 = Task(NCD6, RES_TYPES_LIST, 5, 3)
task5 = Task(NCDA8, RES_TYPES_LIST, 4, 3)
task6 = Task(TSF12, RES_TYPES_LIST, 4, 2)
block_find = BlockFinder([typeX, typeN, typeC], 1, NC2, 1)


def main():
    task6.find_scheme()
    # block_find.find()

if __name__ == "__main__":
    main()
