from classes.constants import PatternClass, Constants
import copy


ITERATION_4_OUTPUT = 10000


class CLOptimizer:

    def __init__(self, parameters, logger):
        self.solution_found = False
        self.logger = logger
        self.samples = parameters["samples"]
        self.optimize_price = parameters["optimize_price"]
        self.iteration = 0
        self.residues2label = []
        self.final_depth = len(self.residues2label)
        self.solutions = 0
        self.solution = Solution(None)
        self.back_up_solutions = []
        self.counter = [0]
        self.depth = 0
        self.symmetry = [self.samples]

    def run(self):
        while not self.solution_found:
            self.find_solution()
            self.samples += 1

    def find_solution(self):
        self.back_up_solutions = []
        self.counter = [0]
        self.depth = 0
        self.generate_residues2label()
        self.create_patterns_codes()
        self.create_solution()
        self.main_cycle()

    def main_cycle(self):
        while True:
            self.next_iteration()

            if self.depth == 0 and self.counter[0] > self.solution.residues[0].patterns_number():
                break

            self.back_up_solutions.append(self.solution.copy())


            if not len(self.residues2label) and self.solution.good:

                if self.optimize_price:
                    pass
                else:
                    pass
                self.solutions += 1

            if self.solution.found and samples_number != self.solution.samples_num:
                break


    def next_iteration(self):
        self.iteration += 1
        # some output

    def go_deeper(self):
        # choose residue with smallest number of patterns
        # add counter
        # depth ++
        pass

    def go_parallel(self):
        # take back-upped scheme
        # pop back-upped scheme
        # counter ++
        pass

    def go_back(self):
        # depth --
        #

    def create_solution(self):
        pass

    def generate_residues2label(self):
        #make a set of residues
        #make residues after and before
        #check good or bad
        #cross out nitrogen patterns
        #if needed calculate prices
        #make labeling array
        #for each residue2label convert patterns to codes
        all_patterns = set().union(*[res.patterns_set for res is self.residues2label])
        all_patterns_list = list(all_patterns)


        self.label_list = sorted(self.label_set, key=lambda word: [alphabet[c] for c in word])


        patterns_codes = PatternsCodes(all_patterns_list)
        self.solution = Solution(patterns_codes)
        pass

    def create_patterns_codes(self):
        pass


class Residue2Label:

    def __init__(self, name, ncs, samples, label_options, residues_after, residues_before):
        self.name = name
        self.samples = samples
        self.pattern_price = dict()
        self.residues_after = residues_after
        self.residues_before = residues_before
        self.patterns_list = []
        self.patterns_codes = []
        self.label_options = label_options

        self.has_15n = False
        for label in label_options:
            if label in Constants.NITRO_TYPES:
                self.has_15n = True
                break

        self.labeling_prices = dict()
        self.patterns_set = self._generate_initial_set(self.samples)
        self._cross_out_N_power(ncs)


    def _generate_initial_set(self, samples):
        # recursive function, that generates all possible combinations
        # of labels given the number of samples for the given residue

        if samples == 0:
            new_set = set()
            new_set.add("")
            return new_set
        current_set = self._generate_initial_set(samples - 1)
        new_set = set()
        for item in current_set:
            for option in self.label_options:
                new_set.add(item + option)
        return new_set

    def _cross_out_N_power(self, ncs):
        if not self.has_15n:
            return
        cross_out_set = set()
        for pattern in self.patterns_set:
            if not self.check_N_power(pattern, ncs):
                cross_out_set.add(pattern)
        self.patterns_set = self.patterns_set.difference(cross_out_set)

    def check_N_power(self, pattern, ncs):
        max_pairs = 1
        got_nitro = False
        for label in pattern:
            max_pairs *= ncs.label_power[label]
            if label in Constants.NITRO_TYPES:
                got_nitro = True
        return got_nitro and max_pairs >= len(self.residues_after)


    def cross_out(self):
        pass

    def patterns_number(self):
        return len(self.patterns_list)

    def cross_out_symmetry(self, symmetry):
        for pattern_code


class PatternsCodes:

    def __init__(self, patterns, ncs):
        self.patterns = patterns
        self.ncs = ncs
        self.codes = []
        self.codes_list = []
        pattern_class = PatternClass()
        self.simple_list = [pattern_class.simplify_pattern(pattern) for pattern in patterns]
        self._create_codes_table()

    def _create_codes_table(self):
        for pattern1 in self.patterns:
            codes_row = []
            for pattern2 in self.patterns:
                symbol_code = self.ncs.calc_code(pattern1, pattern2)
                if symbol_code in self.codes_list:
                    code_number = self.codes_list.index(symbol_code)
                else:
                    code_number = len(self.codes_list)
                    self.codes_list.append(symbol_code)
                codes_row.append(code_number)
            self.codes.append(codes_row)

    def get_pattern_number(self, pattern):
        return self.patterns.index(pattern)

    def get_pattern_by_code(self, pattern_code):
        return self.patterns[pattern_code]

    def get_simple_form_by_code(self, pattern_code):
        return  self.simple_list[pattern_code]

    def get_symbol_code_by_number(self, code_number):
        return self.codes_list[code_number]


class Solution:

    def __init__(self, patterns_codes):
        self.patterns_codes = patterns_codes
        self.patterns = []
        self.residues = []
        self.codes = set()
        self.price = 0
        self.new_codes = set()

    def try_label(self, pattern_code, residue):
        self.new_codes = set()
        if residue.name in residue.residues_after:
            self.new_codes.add(self.patterns_codes.codes[pattern_code][pattern_code])
        for i in range(len(self.residues)):
            res = self.residues[i]
            if res.name in residue.residues_after:
                new_code = self.patterns_codes.codes[pattern_code][self.patterns[i]]
                if new_code in self.codes or new_code in self.new_codes:
                    return False
                else:
                    self.new_codes.add(new_code)
            if res.name in residue.residues_before:
                new_code = self.patterns_codes.codes[self.patterns[i]][pattern_code]
                if new_code in self.codes or new_code in self.new_codes:
                    return False
                else:
                    self.new_codes.add(new_code)
        return True

    def add_label(self, pattern_code, residue):
        self.codes.update(self.new_codes)
        self.new_codes = set()
        self.patterns.append(pattern_code)
        self.residues.append(residue)
        self.price += residue.pattern_price[pattern_code]



