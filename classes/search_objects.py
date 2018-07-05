import sys
import copy
from .constants import Constants


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

    def sort(self):
        for i in range(len(self.patterns)-1):
            for j in range(len(self.patterns)-1-i):
                if self.pattern_bigger(self.patterns[i], self.patterns[i+j+1]):
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

    def pattern_bigger(self, pattern1, pattern2):
        for i in range(len(pattern1)):
            type1 = pattern1[i]
            type2 = pattern2[i]
            if Constants.TYPES.index(type1) > Constants.TYPES.index(type2):
                return True
            if Constants.TYPES.index(type1) < Constants.TYPES.index(type2):
                return False
        return True

    def add_pattern(self, new_pattern):
        if self.try_pattern(new_pattern):
            self.patterns.append(new_pattern)
            self.codes.update(self.new_codes)
            self.simplify()

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
        simple_form = [0 for _ in range(len(Constants.TYPES))]
        for label in pattern:
            for i in range(len(Constants.TYPES)):
                if Constants.TYPES[i] == label:
                    simple_form[i] += 1
                    continue
        result = "".join([str(a) for a in simple_form])
        return result

    def __eq__(self, scheme):
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
        sys.stdout.flush()

    def size(self):
        return len(self.patterns)

    def __str__(self):
        output = ""
        for pattern in self.patterns:
            output +=  pattern + "\n"
        return output

    def __mul__(self, other):
        new_patterns = []
        for pattern_1 in self.patterns:
            for pattern_2 in other.patterns:
                new_patterns.append(pattern_1 + pattern_2)
        new_name = self.name + "_X_" + other.name
        samples = len(new_patterns[0])
        new_scheme = Scheme(new_name, self.ncs, samples, new_patterns)
        return new_scheme


class Product:

    def __init__(self, all_blocks, product_list):
        self.blocks = all_blocks
        self.product_list = product_list
        self.counters = [0 for block_type in self.product_list]
        self.max_counters = [len(self.blocks[block_type[0]][block_type[1]]) for block_type in self.product_list]
        self.overflow = False
        self.ncs = self.blocks[self.product_list[0][0]][self.product_list[0][1]][0].ncs
        self.last_blocks = []

    def __len__(self):
        length = 1
        for counter in self.max_counters:
            length *= counter
        return length

    def __next__(self):
        if self.overflow:
            raise StopIteration
        scheme = Scheme("", self.ncs, 0, [''])
        self.last_blocks = []
        for i in range(len(self.counters)):
            samples_n = self.product_list[i][0]
            patterns_n = self.product_list[i][1]
            curr_block = self.blocks[samples_n][patterns_n][self.counters[i]]
            scheme *= curr_block
            self.last_blocks.append(curr_block)
        for i in range(len(self.counters)):
            self.counters[i] += 1
            if self.counters[i] == self.max_counters[i]:
                self.overflow = True
                self.counters[i] = 0
            else:
                self.overflow = False
                break
        return scheme

    def __iter__(self):
        self.counters = [0 for block_type in self.product_list]
        return self

    def __str__(self):
        return str(self.product_list)


class BestScheme:

    def __init__(self, scheme, blocks, price, label_dict, residues):
        self.scheme = scheme
        self.blocks = blocks
        self.price = price
        self.label_dict = label_dict
        self.samples = scheme.samples
        self.residues = residues
