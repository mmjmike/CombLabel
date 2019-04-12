from classes.constants import PatternClass, Constants
from operator import attrgetter
import copy
import time
import socket


ITERATION_4_OUTPUT = 1000000


class CLOptimizer:

    def __init__(self, parameters, logger, sequence):
        self.solution_found = False
        self.logger = logger
        self.samples = parameters["samples"]
        self.optimize_price = parameters["optimize_price"]
        self.ncs = parameters["ncs"]
        self.sequence = sequence
        self.jobname = parameters["jobname"]
        self.iteration = 0
        self.residues2label = set()
        self.residues_not_label = set()
        self.final_depth = len(self.residues2label)
        self.max_depth = 0
        self.solutions = 0
        self.solution = None
        self.back_up_solutions = []
        self.counter = []
        self.depth = 0
        self.symmetry = []
        self.symmetry.append([self.samples])
        self.best_solution = Solution(None)
        self.patterns_codes = None
        self.current_res = None
        self.t0 = None

    def run(self):
        self.t0 = time.time()

        out = "Search for {} solution started\n".format(self.jobname)
        out += "Host name is:\t\t{}\n".format(socket.gethostname())
        out += "Start clock time: \t\t{}".format(time.strftime("%d-%m-%Y %H:%M:%S", time.gmtime(self.t0)))
        self.logger.info(out)

        while not self.solution_found:
            out = "\nSearch in {} sample(s)".format(self.samples)
            self.logger.info(out)

            self.find_solution()

            out = "Search in {} sample(s) finished".format(self.samples)
            if self.solution_found:
                out += " successfully!"
            else:
                out += ". Max depth reached: {}. Iterations used: {}".format(self.max_depth,
                                                                             self.iteration)
            self.logger.info(out)

            self.samples += 1

        out = "Search for solution finished successfully!\n"
        out += "Calculation time: \t\t{}\n".format(time.strftime("%H:%M:%S",
                                                               time.gmtime(int(time.time() - self.t0))))
        out += "Total iterations: {}".format(self.iteration)
        self.logger.info(out)

    def find_solution(self):
        self.current_res = None
        self.back_up_solutions = []
        self.depth = 0
        self.generate_residues2label()
        self.solution = Solution(self.patterns_codes)
        # choose the first residue to start with
        if not self.residues2label:
            return
        self.symmetry = []
        self.symmetry.append([self.samples])
        self.current_res = min(self.residues2label, key=attrgetter('score'))

        self.current_res.cross_out_symmetry(self.symmetry[-1], self.patterns_codes)
        self.current_res.make_pattern_list()
        self.residues2label.remove(self.current_res)
        self.counter = [0]

        self.main_cycle()

    def main_cycle(self):

        while True:
            # technical update, write data to log and other stuff
            self.next_iteration()

            # if the counter overflows, then we return to the previous depth
            if self.counter[-1] + 1 > self.current_res.patterns_number:
                if self.depth == 0:
                    break
                self.go_back()
                continue

            result, cross_out_task = self.solution.add_label(self.current_res.patterns_list[self.counter[-1]],
                                                             self.current_res)
            # if label is not added, then use next pattern, e.g. go parallel
            if not result:
                self.go_parallel()
                continue

            # if it is the last residue, then we found a solution!
            if not self.residues2label:
                self.solution_found = True

                # in case of price optimization we compare the price of the new solution with
                # the price of the best one
                if self.optimize_price:
                    if not self.solutions:
                        self.solution.copy_solution_to(self.best_solution)
                        self.solutions = 1
                        self.output_solution()
                    elif self.solution.price < self.best_solution.price:
                        self.solution.copy_solution_to(self.best_solution)
                        self.solutions += 1
                        self.output_solution()
                    self.current_res = self.solution.remove_last_label()
                    self.go_parallel()
                    # we don't count solutions without price improvement
                # if we don't optimize price, then after the first solution we break the cycle
                else:
                    self.solutions = 1
                    self.solution.copy_solution_to(self.best_solution)
                    self.output_solution()
                    break
            # if not last residue, we perform the cross out for following residues
            # and go to the next level
            # the case of no patterns left for the following residue is dealt with
            # in the first "if" statement in the cycle
            else:
                self.perform_cross_out(cross_out_task)
                if self.optimize_price and self.solution_found:
                    if not self.check_price():
                        self.restore_last_cross_outs()
                        self.current_res = self.solution.remove_last_label()
                        self.go_parallel()
                        continue
                self.go_deeper()

    def check_price(self):
        price = 0
        for residue in self.residues2label:
            price += residue.find_cheapest_option()
        return price + self.solution.price < self.best_solution.price

    def perform_cross_out(self, cross_out_task):
        for residue in self.residues2label:
            residue.cross_out(cross_out_task, self.patterns_codes)

    def restore_last_cross_outs(self):
        for residue in self.residues2label:
            residue.restore_last_cross_out()

    def next_iteration(self):
        self.iteration += 1

        if self.iteration % ITERATION_4_OUTPUT == 0:

            out = "UPDATE: search in {} samples; iteration {:>9}; time: {:>6d} sec; ".format(self.samples, self.iteration,
                                                             int(time.time() - self.t0))
            out += "max depth={:<2}; solutions found= {:<6}".format(self.max_depth + 1, self.solutions)
            if self.optimize_price and self.solution_found:
                out += "; best price={}".format(self.best_solution.price)
            self.logger.info(out)

    def go_deeper(self):
        self.update_symmetry()
        self.current_res = min(self.residues2label, key=attrgetter('score'))
        self.residues2label.remove(self.current_res)
        self.current_res.cross_out_symmetry(self.symmetry[-1], self.patterns_codes)
        self.current_res.make_pattern_list()
        self.counter.append(0)
        self.depth += 1
        if self.depth > self.max_depth:
            self.max_depth = self.depth
            out = "New max depth ({}) reached at {} samples".format(self.max_depth, self.samples)
            self.logger.info(out)

    def update_symmetry(self):
        current_symmetry = self.symmetry[-1]
        if len(self.symmetry[-1]) == self.samples:
            self.symmetry.append(current_symmetry)
            return

        last_pattern = self.current_res.patterns_list[self.counter[-1]]
        string_pattern = self.patterns_codes.patterns[last_pattern]
        new_symmetry = []
        start_point = 0
        for i in range(len(current_symmetry)):
            if current_symmetry[i] > 1:
                number = 1
                for j in range(current_symmetry[i] - 1):
                    if string_pattern[start_point + j] == string_pattern[start_point + j + 1]:
                        number += 1
                    else:
                        new_symmetry.append(number)
                        number = 1
                    if j == current_symmetry[i] - 2:
                        new_symmetry.append(number)
            else:
                new_symmetry.append(1)
            start_point += current_symmetry[i]
        self.symmetry.append(new_symmetry)

    def go_parallel(self):
        self.counter[self.depth] += 1

    def go_back(self):
        self.depth -= 1
        self.symmetry.pop()
        self.counter.pop()
        self.counter[-1] += 1
        # self.back_up_solutions.pop()
        self.current_res.restore_symmetry()
        self.residues2label.add(self.current_res)
        for residue in self.residues2label:
            residue.restore_last_cross_out()
        self.current_res = self.solution.remove_last_label()

    def generate_residues2label(self):
        self.residues2label = set()
        for res_name in self.sequence.residues:
            residue_obj = self.sequence.residues[res_name]
            if residue_obj.need_label:
                residue_obj.make_patterns(self.samples, self.ncs)
                self.residues2label.add(residue_obj)
            else:
                self.residues_not_label.add(residue_obj)

        all_patterns = set().union(*[res.patterns_set for res in self.residues2label])
        all_patterns_list = list(all_patterns)
        alphabet = Constants.LABEL_SORT_ALPHABET
        patterns = sorted(all_patterns_list, key=lambda word: [alphabet[c] for c in word])

        self.logger.info("Precalculating NMR codes...")
        self.patterns_codes = PatternsCodes(patterns, self.ncs)
        self.logger.info("NMR codes calculated")
        for res in self.residues2label:
            res.translate_patterns(self.patterns_codes)

    def output_solution(self):
        msg = "New solution found! (#{})".format(self.solutions)
        msg += ". Calculation time: {}".format(time.strftime("%H:%M:%S",
                                                        time.gmtime(int(time.time() - self.t0))))
        if self.optimize_price:
            msg += ". New best price = {}".format(str(round(self.best_solution.price, 2)))
        self.logger.info(msg)

        filename = "{}_all_solutions.txt".format(self.jobname)
        out = "[solution]\n"
        out += "% Solution number = {}\n".format(self.solutions)
        out += "% Solution price  = {}\n".format(str(round(self.best_solution.price, 2)))
        out += "Res," + ",".join([" S_" + str(n+1) for n in range(self.samples)])
        out += "\n"
        for i in range(len(self.best_solution.residues)):
            res = self.best_solution.residues[i]
            pattern_number = self.best_solution.patterns[i]
            pattern = self.patterns_codes.get_pattern_by_number(pattern_number)
            out += Constants.TO_THREE_LETTER_CODE[res.name] + ",   "
            out += ",   ".join([char for char in pattern])
            out += "\n"
        non_labeled_residues = []

        for res in self.sequence.residues:
            res_obj = self.sequence.residues[res]
            if not res_obj.need_label:
                three_letter_name = Constants.TO_THREE_LETTER_CODE[res]
                non_labeled_residues.append(three_letter_name)

        for res_name in non_labeled_residues:
            out += res_name + ",   X" * self.samples
        out += "\n\n"
        mode = "w"
        if self.solutions > 1:
            mode = "a"
        with open(filename, mode=mode) as f:
            f.write(out)
        self.logger.debug(out)
        self.write_dict_and_stats()

    def write_dict_and_stats(self):
        solution = self.best_solution.dict_form
        solution = add_unlabeled_residues(solution, self.sequence)
        output_solution = self.jobname + "_solution.txt"
        output_dictionary = self.jobname + "_dictionary.txt"
        label_options = get_stock_from_solution(solution)
        pairs_table = PairsTable(self.sequence, label_options)

        price_dict = self.sequence.get_prices()


        params = {
            "ncs": self.ncs,
            "name": self.jobname,
            "label_options": label_options,
            "sequence": self.sequence,
            "solution": solution,
            "pairs_table": pairs_table,
            "show_price": self.optimize_price,
            "price_dict": price_dict,
            "output_solution": output_solution,
            "output_dictionary": output_dictionary
        }
        write_solution_stats(params, cl=True)


class Residue2Label:

    def __init__(self, name, label_options, prices):
        self.name = name
        self.samples = 0
        self.pattern_price = dict()
        self.patterns_list = []
        self.label_options = label_options
        self.symmetry_cross_out = set()
        self.crossed_out_sets = []
        self.include_other = True
        self.residues_after = set()
        self.residues_before = set()
        self.pattern_codes_set = set()
        self.has_15n = False
        for label in label_options:
            if label in Constants.NITRO_TYPE_NAMES:
                self.has_15n = True
                break
        if not self.has_15n:
            self.include_other = False
        self.has_13c = False
        for label in label_options:
            if label in Constants.CARBON_TYPE_NAMES:
                self.has_13c = True
                break
        self.need_label = False
        if self.has_15n or self.has_13c:
            self.need_label = True
        self.has_double_label = False
        for label in label_options:
            if label in Constants.CARBON_TYPE_NAMES and Constants.NITRO_TYPE_NAMES:
                self.has_double_label = True
                break
        self.labeling_prices = prices
        if not prices:
            for label in label_options:
                self.labeling_prices[label] = 0
        self.patterns_set = set()

    def translate_patterns(self, patterns_codes):
        self.pattern_codes_set = set()
        self.pattern_price = dict()
        self.crossed_out_sets = []
        for pattern in self.patterns_set:
            pattern_code = patterns_codes.patterns.index(pattern)
            price = self.calculate_pattern_price(pattern)
            self.pattern_codes_set.add(pattern_code)
            self.pattern_price[pattern_code] = price

    def calculate_pattern_price(self, pattern):
        price = 0
        for label in pattern:
            price += self.labeling_prices[label]
        return price

    def make_patterns(self, samples, ncs):
        self.samples = samples
        self.patterns_set = self.generate_initial_set(self.samples)
        self.cross_out_N_power(ncs)

    def generate_initial_set(self, samples):
        # recursive function, that generates all possible combinations
        # of labels given the number of samples for the given residue

        if samples == 0:
            new_set = set()
            new_set.add("")
            return new_set
        current_set = self.generate_initial_set(samples - 1)
        new_set = set()
        for item in current_set:
            for option in self.label_options:
                new_set.add(item + option)
        return new_set

    def add_residue_after(self, res_name):
        if self.has_13c:
            self.residues_after.add(res_name)

    def add_residue_before(self, res_name):
        if self.has_15n:
            self.residues_before.add(res_name)

    def cross_out_N_power(self, ncs):
        if not self.has_15n:
            return
        if "Other" in self.residues_before and not self.include_other:
            self.residues_before.remove("Other")
        cross_out_set = set()
        for pattern in self.patterns_set:
            if not self.check_N_power(pattern, ncs):
                cross_out_set.add(pattern)
        self.patterns_set = self.patterns_set.difference(cross_out_set)

    def check_N_power(self, pattern, ncs):
        max_pairs = 1
        got_nitro = False
        for label in pattern:
            label_type = Constants.BASIC_TYPES[Constants.TYPES_NAMES.index(label)]
            label_power = ncs.label_power[label_type]
            max_pairs *= label_power
            if label in Constants.NITRO_TYPE_NAMES:
                got_nitro = True
        return got_nitro and max_pairs >= len(self.residues_before)

    def cross_out(self, task, patterns_codes):
        cross_out_set = set()
        for residue in task:
            if residue in self.residues_after:
                pattern_after = task[residue][0]
                codes = set(task[residue][1:])
                for pattern in self.pattern_codes_set:
                    if patterns_codes.codes[pattern][pattern_after] in codes:
                        cross_out_set.add(pattern)
        self.crossed_out_sets.append(cross_out_set)
        self.pattern_codes_set = self.pattern_codes_set.difference(cross_out_set)

    @property
    def patterns_number(self):
        return len(self.pattern_codes_set)

    @property
    def score(self):
        return len(self.pattern_codes_set)*1000 - ord(self.name)

    def cross_out_symmetry(self, symmetry, pattern_codes):
        self.symmetry_cross_out = set()
        if len(symmetry) == self.samples:
            return
        label_num_dict = Constants.LABEL_SORT_ALPHABET
        for pattern_code in self.pattern_codes_set:
            pattern = pattern_codes.get_pattern_by_number(pattern_code)
            start_point = 0
            for i in range(len(symmetry)):
                if symmetry[i] > 1:
                    for j in range(symmetry[i] - 1):
                        if (label_num_dict[pattern[start_point + j]]
                                > label_num_dict[pattern[start_point + j + 1]]):
                            self.symmetry_cross_out.add(pattern_code)
                start_point += symmetry[i]
        self.pattern_codes_set = self.pattern_codes_set.difference(self.symmetry_cross_out)

    def restore_symmetry(self):
        self.pattern_codes_set = self.pattern_codes_set.union(self.symmetry_cross_out)

    def restore_last_cross_out(self):
        last_cross_out = self.crossed_out_sets.pop()
        self.pattern_codes_set = self.pattern_codes_set.union(last_cross_out)

    def make_pattern_list(self):
        self.patterns_list = []
        self.patterns_list = sorted(self.pattern_codes_set)

    def find_cheapest_option(self):
        if not self.pattern_codes_set:
            return 0
        else:
            patterns_prices = [self.pattern_price[pattern_code] for pattern_code in self.pattern_codes_set]
            return min(patterns_prices)


class PatternsCodes:

    def __init__(self, patterns, ncs):
        self.patterns = patterns
        if patterns:
            self.samples = len(patterns[0])
        else:
            self.samples = 0
        other_pattern = "X"*self.samples
        if other_pattern not in self.patterns:
            self.patterns.append(other_pattern)
        self.ncs = ncs
        self.codes = []
        self.codes_list = []
        pattern_class = PatternClass()
        self.simple_list = [pattern_class.simplify_pattern(pattern) for pattern in patterns]
        self._create_codes_table()
        self.other_code = self.patterns.index(other_pattern)

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

    def get_pattern_by_number(self, pattern_number):
        return self.patterns[pattern_number]

    def get_simple_form_by_number(self, pattern_number):
        return  self.simple_list[pattern_number]

    def get_symbol_code_by_number(self, code_number):
        return self.codes_list[code_number]


class Solution:

    def __init__(self, patterns_codes):
        self.patterns_codes = patterns_codes
        self.patterns = []
        self.residues = []
        self.codes = set()
        self.price = 0
        self.found = False
        self.new_codes = set()
        self.back_up_codes = []
        self.price_back_up = []

    def add_label(self, pattern_code, residue):
        # task: which res, which pattern second, code value
        cross_out_task = dict()
        self.new_codes = set()
        residue_added = False

        if residue.name in residue.residues_after:
            new_code = self.patterns_codes.codes[pattern_code][pattern_code]
            if new_code in self.codes:
                return False, cross_out_task
            self.new_codes.add(new_code)
            cross_out_task.update({residue.name: [pattern_code, new_code]})
            residue_added = True
        if "Other" in residue.residues_before and residue.include_other:
            other_code = self.patterns_codes.other_code
            new_code = self.patterns_codes.codes[other_code][pattern_code]
            if new_code in self.codes or new_code in self.new_codes:
                return False, cross_out_task
            if not residue_added:
                cross_out_task.update({residue.name: [pattern_code, new_code]})
                residue_added = True
            else:
                cross_out_task[residue.name].append(new_code)
            self.new_codes.add(new_code)
        for i in range(len(self.residues)):
            res = self.residues[i]
            if res.name in residue.residues_after:
                new_code = self.patterns_codes.codes[pattern_code][self.patterns[i]]
                if new_code in self.codes or new_code in self.new_codes:
                    return False, cross_out_task
                self.new_codes.add(new_code)
                cross_out_task.update({res.name: [self.patterns[i], new_code]})
            if res.name in residue.residues_before:
                new_code = self.patterns_codes.codes[self.patterns[i]][pattern_code]
                if new_code in self.codes or new_code in self.new_codes:
                    return False, cross_out_task
                self.new_codes.add(new_code)
                if not residue_added:
                    cross_out_task.update({residue.name: [pattern_code, new_code]})
                else:
                    cross_out_task[residue.name].append(new_code)

        self.codes.update(self.new_codes)
        self.back_up_codes.append(self.new_codes)
        self.new_codes = set()
        self.patterns.append(copy.copy(pattern_code))
        self.residues.append(residue)
        self.price_back_up.append(copy.copy(self.price))
        self.price += residue.pattern_price[pattern_code]
        return True, cross_out_task

    def remove_last_label(self):
        self.patterns.pop()
        self.price = self.price_back_up.pop()
        self.codes = self.codes.difference(self.back_up_codes.pop())
        return self.residues.pop()

    def copy_solution_to(self, another):
        another.patterns = copy.copy(self.patterns)
        another.residues = copy.copy(self.residues)
        another.price = copy.copy(self.price)
        another.found = copy.copy(self.found)
        another.patterns_codes = copy.copy(self.patterns_codes)
        another.new_codes = copy.copy(self.new_codes)
        another.codes = copy.copy(self.codes)

    @property
    def dict_form(self):
        solution_dict = {}
        for i in range(len(self.residues)):
            res = self.residues[i].name
            pattern = self.patterns_codes.get_pattern_by_number(self.patterns[i])
            solution_dict[res] = pattern
        return solution_dict

    def __str__(self):
        out = ""
        for i in range(len(self.residues)):
            res_name = self.residues[i].name
            pattern = self.patterns_codes.get_pattern_by_number(self.patterns[i])
            out += "{}: {}\n".format(res_name, pattern)
        return out


class PairsTable:

    def __init__(self, sequence, label_options):
        self.label_options = label_options
        self.sequence = sequence.sequence
        self.unique_pairs = []
        self.residues_first = []
        self.residues_second = []
        self.res_no_first = list(Constants.RES_TYPES_LIST)
        self.res_no_second = list(Constants.RES_TYPES_LIST)
        self.residue_pairs = []
        self.all_residue_pairs = []
        self.bad_residues = []
        self.min_nitrogens = []
        self.residues_nitro = []  # list of residues with 15N labels. Last one can be "Other"
        self.residues_not_nitro = []  # list of residues without 15N labels
        self.residues_carbon = []  # list of residues with 13C labels. Last one can be "Other"
        self.residues_not_carbon = []  # list of residues without 13C labels
        self.residues_to_label = []  # ordered list of residues with 13C or 15N labels
        self.non_labeled_residues = []  # residues without any labels
        self.res_has_diagonal = {}
        self._calculate_stats()

    def _calculate_stats(self):
        # Stock class object required for calculating stats
        # Almost all methods should be rewritten as some of class variables
        # are defined and created in the methods other than '__init__'

        self._rank_residues()       # rank all residues in sequence by the number of
                                    # distinct pairs (two or more equal pairs are counted as one)
        self._residues_to_label()   # creating a list of residues that are to be labeled based
                                    # on the residue ranks and presence of labels in stock
        self._count_pairs()         # counting all the pairs present in sequence
        self._check_carbon_pairs()  # recount pairs in order to skip empty rows
        self._count_nitrogens()     # counting non-epty cells in columns
                                    # for residues, that have 15N-labelin in stock
        self._find_bad_residues()   # find residues, that have diagonal pair and pair with
                                    # "Other" residues and don't have C/N labels in stock
                                    # so that you can't tell diagonal cell from pair with "Other"
        self._calculate_subtable_coordinates()      # precalculate coordinates of
                                                    # subtable cells for each residue
        self._calculate_under_cells()   # find cells that are under the particular one for crossing-out

    def _rank_residues(self):
        # Bad code for initializing outside __init__
        residue_rank = [0 for _ in Constants.RES_TYPES_LIST]
        self.unique_pairs = []
        for i in range(len(self.sequence) - 1):
            pair = self.sequence[i:i+2]
            if pair[0] in Constants.RES_TYPES_LIST and pair[1] in Constants.RES_TYPES_LIST:
                if pair not in self.unique_pairs:
                    self.unique_pairs.append(pair)
                    if pair[1] != pair[0]:
                        residue_rank[Constants.RES_TYPES_LIST.index(pair[1])] += 1
                    residue_rank[Constants.RES_TYPES_LIST.index(pair[0])] += 1
                if pair[0] in self.res_no_first:
                    self.res_no_first.pop(self.res_no_first.index(pair[0]))
                if pair[1] in self.res_no_second:
                    self.res_no_second.pop(self.res_no_second.index(pair[1]))
        self.ranked_residues, self.rank_of_residue = sort_residues(Constants.RES_TYPES_LIST, residue_rank)
        self.res_no_second.append("P")

        # residues that are first in some pair in sequence
        self.residues_first = [res for res in self.ranked_residues if res not in self.res_no_first]

        # residues that are second in some pair in sequence
        self.residues_second = [res for res in self.ranked_residues if res not in self.res_no_second]

    def _count_pairs(self):
        # table of residue pairs taken in account labeling stock
        self.residue_pairs = [[0 for col in range(len(self.residues_nitro))]
                              for row in range(len(self.residues_carbon))]
        # table of all residue pairs present in protein sequence
        self.all_residue_pairs = [[0 for col in range(len(self.residues_second))]
                                  for row in range(len(self.residues_first))]
        for i in range(len(self.sequence) - 1):  # Count all pairs in sequence,
            pair = self.sequence[i:(i + 2)]
            if pair[0] == '*' or pair[1] == '*':
                continue
            if pair[0] not in self.residues_carbon:
                first_res = len(self.residues_carbon) - 1
            else:
                first_res = self.residues_carbon.index(pair[0])
            if pair[1] not in self.residues_nitro:
                second_res = len(self.residues_nitro) - 1
            else:
                second_res = self.residues_nitro.index(pair[1])
            self.residue_pairs[first_res][second_res] += 1
            if pair[0] == pair[1]:
                self.res_has_diagonal[pair[0]] = True
            index_first = self.residues_first.index(pair[0])
            if pair[1] != "P":
                index_second = self.residues_second.index(pair[1])
                self.all_residue_pairs[index_first][index_second] += 1

    def _residues_to_label(self):
        for i in range(len(self.ranked_residues)):
            residue = self.ranked_residues[i]
            if self.rank_of_residue[i]:
                nitro = False
                for type in Constants.NITRO_TYPES:
                    if type in self.label_options[residue]:
                        nitro = True
                        break
                if nitro and residue in self.residues_second:
                    self.residues_nitro.append(residue)
                else:
                    self.residues_not_nitro.append(residue)
                carbon = False
                for type in Constants.CARBON_TYPES:
                    if type in self.label_options[residue]:
                        carbon = True
                        break
                if carbon and residue in self.residues_first:
                    self.residues_carbon.append(residue)
                else:
                    self.residues_not_carbon.append(residue)
                if residue in self.residues_carbon or residue in self.residues_nitro:
                    self.residues_to_label.append(residue)
                    self.res_has_diagonal[residue] = False
                else:
                    self.non_labeled_residues.append(residue)
        if len(self.residues_not_nitro):        # add "Other" column if there
            self.residues_nitro.append("Other") # are residues not labeled by 15N
        if len(self.residues_not_carbon):       # add "Other" row if there
            self.residues_carbon.append("Other")# are residues not labeled by 13C

    def _calculate_subtable_coordinates(self):
        # precalculate coordinates of
        # subtable cells for each residue

        self.subtable_coordinates = []
        check_other = (self.residues_carbon and self.residues_carbon[-1] == 'Other')
        carbon = 0
        nitro = 0
        total_cells = 0

        for i in range(len(self.residues_to_label)):
            self.new_coordinates = []
            residue = self.residues_to_label[i]
            if i == 0:
                if residue in self.residues_nitro:
                    if residue in self.residues_carbon:
                        if self.residue_pairs[0][0]:
                            self.new_coordinates.append((0, 0))
                        carbon += 1
                    if check_other and self.residue_pairs[-1][0] and residue not in self.bad_residues:
                        self.new_coordinates.append((len(self.residues_carbon) - 1, 0))
                    nitro += 1
            else:
                if residue in self.residues_carbon:
                    for j in range(nitro):
                        if self.residue_pairs[carbon][j]:
                            self.new_coordinates.append((carbon, j))
                    carbon += 1
                if residue in self.residues_nitro:
                    for j in range(carbon):
                        if self.residue_pairs[j][nitro]:
                            self.new_coordinates.append((j, nitro))
                    if check_other and self.residue_pairs[-1][nitro] and residue not in self.bad_residues:
                        self.new_coordinates.append((len(self.residues_carbon) - 1, nitro))
                    nitro += 1
            total_cells += len(self.new_coordinates)
            self.subtable_coordinates.append(self.new_coordinates)
        curr_cells = 0
        for i in range(len(self.subtable_coordinates)):
            add = len(self.subtable_coordinates[i])
            curr_cells += add

    def _calculate_under_cells(self):
        # find cells that are under the particular one for crossing-out

        self.under_cells = {}
        got_other = 0
        if self.residues_carbon and self.residues_carbon[-1] == "Other":
            got_other = 1
        rows = len(self.residues_carbon) - got_other
        for cell_list in self.subtable_coordinates:
            for cell in cell_list:
                under_cells = []
                row = cell[0] + 1
                if row > rows:
                    row = 0     # if the cell is in "Other" row,
                                #  then all other cells in the column are under it
                column = cell[1]
                for i in range(rows - row):
                    if self.residue_pairs[row+i][column]:
                        under_cells.append((row+i, column))
                self.under_cells.update({cell: under_cells})

    def _check_carbon_pairs(self):
        carb_other = 0
        changed = False
        if self.residues_carbon and self.residues_carbon[-1] == "Other":
            carb_other = 1
        if self.residues_nitro and self.residues_nitro[-1] == "Other":
            carbons = len(self.residues_carbon) - carb_other
            for i in range(carbons):
                k = carbons - 1 - i
                got_pair = False
                for j in range(len(self.residues_nitro) - 1):
                    if self.residue_pairs[k][j]:
                        got_pair = True
                        break
                if not got_pair:
                    self.residues_not_carbon.append(self.residues_carbon[k])
                    residue = self.residues_carbon[k]
                    del self.residues_carbon[k]
                    if residue not in self.residues_nitro:
                        del self.residues_to_label[self.residues_to_label.index(residue)]
                    changed = True
        self.ranks = {}
        for res in self.residues_to_label:
            self.ranks.update({res: 0})
        if changed:
            self._count_pairs()
        for i in range(len(self.residues_nitro)):
            if self.residues_nitro[i] == "Other":
                break
            for j in range(len(self.residues_carbon)):
                if self.residue_pairs[j][i]:
                    if self.residues_carbon[j] != "Other":
                        self.ranks[self.residues_carbon[j]] += 1
                    if self.residues_carbon[j] != self.residues_nitro[i]:
                        self.ranks[self.residues_nitro[i]] += 1
        ranks = []
        residues = []
        for residue in self.residues_to_label:
            if residue in self.ranks:
                residues.append(residue)
                ranks.append(self.ranks[residue])
        self.residues_to_label, self.residue_ranks = sort_residues(residues, ranks)
        self.residues_nitro_new = []
        self.residues_carbon_new = []
        for res in self.residues_to_label:
            if res in self.residues_carbon:
                self.residues_carbon_new.append(res)
            if res in self.residues_nitro:
                self.residues_nitro_new.append(res)
        if "Other" in self.residues_carbon:
            self.residues_carbon_new.append("Other")
        if "Other" in self.residues_nitro:
            self.residues_nitro_new.append("Other")
        self.residues_carbon = self.residues_carbon_new
        self.residues_nitro = self.residues_nitro_new
        self._count_pairs()

    def _add_rank(self, residue):
        self.ranks[self.residues_to_label.index(residue)] += 1

    def _count_nitrogens(self):
        # counting non-empty cells in columns
        # for residues, that have 15N-labeling in stock
        for i in range(len(self.residues_nitro)):
            pairs_count = 0
            for j in range(len(self.residues_carbon)):
                if self.residue_pairs[j][i]:
                    pairs_count += 1
            self.min_nitrogens.append(pairs_count)

    def _find_bad_residues(self):
        # find residues, that have diagonal pair and pair with
        # "Other" residues and don't have C/N labels in stock
        # so that you can't tell diagonal cell code from code of pair with "Other"

        for i in range(len(self.residues_nitro)):
            res = self.residues_nitro[i]
            if res == "Other":
                continue
            if ("D" not in self.label_options[res]
                  and "S" not in self.label_options[res]
                  and "T" not in self.label_options[res]
                  and self.residues_carbon[-1] == 'Other'
                  and self.residues_nitro[i] in self.residues_carbon
                  and self.residue_pairs[self.residues_carbon.index(res)][i]
                  and self.residue_pairs[-1][i]):
                self.bad_residues.append(res)

    def make_full_pairs_table(self):
        output  = "#"*50+"\n"
        output += "# Table of all amino acid pairs \n"
        output += "# \n"
        output += "# Number in the table represents how many times \n"
        output += "# the pair occurs in the sequence \n\n"
        output += "[full_pairs_table]\n"
        output += "Res," + ",".join([Constants.TO_THREE_LETTER_CODE[res] for res in self.residues_second]) + "\n"
        for i in range(len(self.residues_first)):
            output += Constants.TO_THREE_LETTER_CODE[self.residues_first[i]]
            for j in range(len(self.residues_second)):
                output += "," + "{:>3}".format(self.all_residue_pairs[i][j])
            if i+1 < len(self.residues_first):
                output += "\n"
        return output

    def make_stock_pairs_table(self):
        output = "\n\n" + "#" * 50 + "\n"
        output += "# Table of amino acid pairs used for \n"
        output += "# combinatorial labeling \n"
        output += "# \n"
        output += "# Number in the table represents how many times \n"
        output += "# the pair occurs in the sequence \n\n"
        output += "[pairs_table]\n"
        additional_output = "\n"
        if self.residues_carbon[-1] == "Other":
            output += "   "
        output += "Res,"
        if self.residues_nitro[-1] == "Other":
            output += ",".join([Constants.TO_THREE_LETTER_CODE[res] for res in self.residues_nitro[:-1]])
            output += ",OtherN"
            additional_output += "\nOtherN: " + ",".join(
                [Constants.TO_THREE_LETTER_CODE[res] for res in self.residues_not_nitro])
        else:
            output += ",".join([Constants.TO_THREE_LETTER_CODE[res] for res in self.residues_nitro])
        output += "\n"
        for i in range(len(self.residues_carbon)):
            res1 = self.residues_carbon[i]
            if res1 == "Other":
                output += "OtherC"
                additional_output += "\nOtherC: " + ",".join(
                    [Constants.TO_THREE_LETTER_CODE[res] for res in self.residues_not_carbon])
            else:
                if self.residues_carbon[-1] == "Other":
                    output += "   "
                output += Constants.TO_THREE_LETTER_CODE[res1]
            for j in range(len(self.residues_nitro)):

                output += "," + "{:>3}".format(self.residue_pairs[i][j])
            if i + 1 < len(self.residues_carbon):
                output += "\n"
        output += additional_output + "\n"
        return output

    def make_pairs_codes(self, solution, ncs):
        samples_num = len(next(iter(solution.values())))

        output = "\n\n" + "#" * 50 + "\n"
        output += "# Spectrum codes of the labeled amino acid pairs \n#\n"
        output += "# Amino acid and labeling code strings \n"
        output += "# according to the number of samples are in the headers\n"
        output += "# Spectrum codes are in the table\n\n"
        output += "[pairs_codes]\n"

        additional_output = "\n"
        separator = ","
        if samples_num > 3:
            separator += " " * (samples_num - 3)
        otherC_spaces = ""
        if self.residues_carbon[-1] == "Other":
            otherC_spaces += "   "
        output += otherC_spaces + "   ," + " " * samples_num + separator
        if self.residues_nitro[-1] == "Other":
            res_list = self.residues_nitro[:-1]
        else:
            res_list = self.residues_nitro
        output += separator.join([Constants.TO_THREE_LETTER_CODE[res] for res in res_list])
        if self.residues_nitro[-1] == "Other":
            output += ",OtherN"
            additional_output += "\nOtherN: " + ",".join(
                [Constants.TO_THREE_LETTER_CODE[res] for res in self.residues_not_nitro])
        output += "\n"
        output += otherC_spaces + "   ," + " " * samples_num
        for res in res_list:
            output += "," + solution[res]
        if self.residues_nitro[-1] == "Other":
            output += "," + "X" * samples_num
        output += "\n"
        for i in range(len(self.residues_carbon)):
            res1 = self.residues_carbon[i]
            if res1 == "Other":
                additional_output += "\nOtherC: " + ",".join(
                    [Constants.TO_THREE_LETTER_CODE[res] for res in self.residues_not_carbon])
                output += "OtherC"
            else:
                output += otherC_spaces + Constants.TO_THREE_LETTER_CODE[res1]
            output += "," + solution[res1]
            for j in range(len(self.residues_nitro)):
                res2 = self.residues_nitro[j]
                if res2 == "Other":
                    code = "0" * samples_num
                elif self.residue_pairs[i][j]:
                    code = ncs.calc_code(solution[res1], solution[res2])
                else:
                    code = " " * samples_num
                output += "," + code
            if i + 1 < len(self.residues_carbon):
                output += "\n"
        output += additional_output
        output += "\n"
        return output


def calculate_code_dictionary(solution, sequence, ncs):
    codes_table = []
    meanings_table = []
    code_meanings = {}

    for i in range(len(sequence.sequence) - 1):
        pair = sequence.sequence[i:i+2]
        residue1 = pair[0]
        residue2 = pair[1]
        pattern1 = solution[residue1]
        pattern2 = solution[residue2]
        code = ncs.calc_code(pattern1, pattern2)
        meaning = "".join([residue1, str(i+1), " - ", residue2, str(i+2)])
        codes_table.append(code)
        meanings_table.append(meaning)
        if code not in code_meanings:
            code_meanings[code] = [pair]
        else:
            code_meanings[code].append(pair)

    sorted_meanings = [mean for code, mean in sorted(zip(codes_table,
                                                            meanings_table))]
    sorted_codes_table = [code for code in sorted(codes_table)]
    return sorted_codes_table, sorted_meanings, code_meanings


def add_unlabeled_residues(solution, sequence):
    res_types = sequence.residue_types_used
    samples = len(next(iter(solution.values())))
    for residue in res_types:
        if residue not in solution:
            solution[residue] = "X" * samples
    solution["Other"] = "X" * samples
    return solution


def write_solution_stats(parameters, cl=False):
    sequence = parameters["sequence"]
    solution = parameters["solution"]
    ncs = parameters["ncs"]

    codes, meanings, meanings_dict = calculate_code_dictionary(solution, sequence, ncs)
    write_code_dict(parameters, codes, meanings, cl=cl)

    stats = calculate_stats(parameters, meanings_dict)
    pairs_table = parameters["pairs_table"]

    output = ""
    output += pairs_table.make_full_pairs_table()
    output += pairs_table.make_stock_pairs_table()
    output += "\n\n"
    output += str(ncs)
    output += make_solution_output(parameters)
    output += pairs_table.make_pairs_codes(solution, ncs)
    output += stats_to_text(stats)
    with open(parameters["output_solution"], "w") as f:
        f.write(output)
        f.close()
    if not cl:
        print("Solution stats written to file '{}'".format(parameters["output_solution"]))


def write_code_dict(parameters, codes, meanings, cl=False):
    output = "#" * 50 + "\n"
    output += "# The Spectrum codes dictionary for \'" + parameters["name"] + "\':\n"
    output += "# Spectrum code: First AA - Second AA\n\n"
    for i in range(len(codes)):
        output += "{}: {}\n".format(codes[i], meanings[i])
    with open(parameters["output_dictionary"], "w") as f:
        f.write(output)
        f.close()
    if not cl:
        print("Code dictionary written to file '{}'".format(parameters["output_dictionary"]))


def make_solution_output(parameters):
    solution = parameters["solution"]
    pairs_table = parameters["pairs_table"]
    output = "[solution]\n"
    output += "% Soluiton number = 0\n"
    if parameters["show_price"]:
        price = 0
        for res in solution:
            if res == "Other":
                continue
            price += calc_price(parameters["price_dict"], res, solution[res])
        output += "% Solution price  = {}\n".format(str(round(price, 2)))
    output += "Res"
    samples_num = len(next(iter(solution.values())))
    for i in range(samples_num):
        output += ",S" + str(i + 1)
    output += "\n"
    for residue in pairs_table.residues_to_label:
        output += Constants.TO_THREE_LETTER_CODE[residue] + ", " + ", ".join(
            list(solution[residue])) + "\n"
    for res in pairs_table.non_labeled_residues:
        output += Constants.TO_THREE_LETTER_CODE[res] + ", " + ", ".join(list("X" * samples_num)) + "\n"
    output += "\n"
    return output


def calculate_stats(parameters, meanings_dict):
    stats = {}
    sequence = parameters["sequence"]
    pairs_table = parameters["pairs_table"]

    N_aa = len(sequence.sequence)
    N_pairs = len(pairs_table.unique_pairs)
    PI_p = len([pair for pair in pairs_table.unique_pairs if pair[1] == "P"])
    PI_a = len([res for res in sequence.sequence if res == 'P'])
    PU_p = 0
    PN_a = 0
    PN_p = 0
    for i in range(len(pairs_table.residues_first)):
        for j in range(len(pairs_table.residues_second)):
            number_of_pairs = pairs_table.all_residue_pairs[i][j]
            if number_of_pairs == 1:
                PU_p += 1
            elif number_of_pairs > 1:
                PN_p += 1
                PN_a += number_of_pairs
    PU_a = PU_p
    SI_a = 0
    SI_p = 0
    SU_p = 0
    SN2_a = 0
    SN2_p = 0
    SN1_a = 0
    SN1_p = 0
    for i in range(len(pairs_table.residues_carbon)):
        for j in range(len(pairs_table.residues_nitro)):
            cell_value = pairs_table.residue_pairs[i][j]
            res1 = pairs_table.residues_carbon[i]
            res2 = pairs_table.residues_nitro[j]
            if cell_value:
                seq = sequence.sequence
                if res2 == 'Other':
                    diff_pairs = len(set([seq[s + 1] for s in range(len(seq) - 1)
                                          if seq[s] == res1
                                          and seq[s + 1] in pairs_table.residues_not_nitro]))
                    if res1 == 'Other':
                        diff_pairs = len(set([seq[s + 1] for s in range(len(seq) - 1)
                                              if seq[s] in pairs_table.residues_not_carbon
                                              and seq[s + 1] in pairs_table.residues_not_nitro]))
                    SI_a += cell_value
                    SI_p += diff_pairs
                elif res1 == 'Other':
                    diff_pairs = len(set([seq[s] for s in range(len(seq) - 1)
                                          if seq[s + 1] == res2
                                          and seq[s] in pairs_table.residues_not_carbon]))
                    if res2 in pairs_table.bad_residues:
                        SN1_p += diff_pairs
                        SN1_a += cell_value
                    elif cell_value == 1:
                        SU_p += 1
                    elif diff_pairs > 1:
                        SN1_p += diff_pairs
                        SN1_a += cell_value
                    else:
                        SN2_p += 1
                        SN2_a += cell_value
                elif (res1 == res2
                      and res2 in pairs_table.bad_residues):
                    SN1_p += 1
                    SN1_a += cell_value
                elif cell_value == 1:
                    SU_p += 1
                else:
                    SN2_p += 1
                    SN2_a += cell_value
    SU_a = SU_p
    stats.update({
        "SIa": SI_a,
        "SIp": SI_p,
        "SUa": SU_a,
        "SUp": SU_p,
        "SN1a": SN1_a,
        "SN1p": SN1_p,
        "SN2a": SN2_a,
        "SN2p": SN2_p,
        "Na": N_aa,
        "Np": N_pairs,
        "PIa": PI_a,
        "PIp": PI_p,
        "PUa": PU_a,
        "PUp": PU_p,
        "PNa": PN_a,
        "PNp": PN_p
    })

    LI_a = 0
    LI_p = 0
    LN1_a = 0
    LN1_p = 0
    LN2_a = 0
    LN2_p = 0
    LU_a = 0
    LA_p = 0
    LA_a = 0

    for code in meanings_dict:
        meaning = meanings_dict[code]
        if set(list(code)) == set("0"):
            LI_p += len(set(meaning))
            LI_a += len(meaning)
        elif len(meaning) == 1:
            LU_a += 1
        elif len(set(meaning)) == 1:
            LN2_p += 1
            LN2_a += len(meaning)
        elif len(set([pair[1] for pair in meaning])) == 1:  # all second residues in pair are equal
            LN1_p += len(set(meaning))
            LN1_a += len(meaning)
        else:
            LA_a += len(meaning)
            LA_p += len(set(meaning))
    LU_p = LU_a

    stats.update({
        "LIa": LI_a,
        "LIp": LI_p,
        "LUa": LU_a,
        "LUp": LU_p,
        "LN1a": LN1_a,
        "LN1p": LN1_p,
        "LN2a": LN2_a,
        "LN2p": LN2_p,
        "LAa": LA_a,
        "LAp": LA_p
    })
    return stats


def stats_to_text(stats):
    output = "\n\n" + "#" * 50 + "\n"
    output += "# Calculation statistics\n"
    output += "# \n\n"
    output += "[stats]\n"
    output += ""
    output += "\n\n"
    output += "# Statistics for PAIRS in amino acid sequence\n"
    output += "# The STOCK (availability of isotopically labeled amino acid)\n"
    output += "# is NOT accounted for in this statistics\n"
    output += "# The labeling scheme is NOT accounted too\n"
    output += "\n"
    output += "[stats,pairs]\n"
    output += "Par, Description,  Residues, Pairs\n"
    output += "N,   Total number, {:>8}, {:>5}\n".format(stats["Na"], stats["Np"])
    output += "PI,  Invisible,    {:>8}, {:>5}\n".format(stats["PIa"], stats["PIp"])
    output += "PU,  Unique,       {:>8}, {:>5}\n".format(stats["PUa"], stats["PUp"])
    output += "PN,  Non-unique,   {:>8}, {:>5}\n".format(stats["PNa"], stats["PNp"])
    output += "\n"
    output += "# Statistics for STOCK-available pairs in amino acid sequence\n"
    output += "# The STOCK is used to check whether the particular pairs are distinguishable \n"
    output += "# in principle with some labeling scheme unlimited in size with some NMR spectra\n"
    output += "# The particular labeling scheme, found by the program, is NOT accounted here\n"
    output += "\n"
    output += "[stats,stock]\n"
    output += "Par, Description,     Residues, Pairs\n"
    output += "N,   Total number,    {:>8}, {:>5}\n".format(stats["Na"], stats["Np"])
    output += "SI,  Invisible,       {:>8}, {:>5}\n".format(stats["SIa"], stats["SIp"])
    output += "SU,  Unique code,     {:>8}, {:>5}\n".format(stats["SUa"], stats["SUp"])
    output += "SN2, AA type of both, {:>8}, {:>5}\n".format(stats["SN2a"], stats["SN2p"])
    output += "SN1, AA type of last, {:>8}, {:>5}\n".format(stats["SN1a"], stats["SN1p"])
    output += "\n"
    output += "# Statistics for LABELING CODES\n"
    output += "# The pairs are distinguishable, if their labeling codes are different\n"
    output += "# Both sequence, stock, NMR spectra and particular labeling scheme is accounted here\n"
    output += "\n"
    output += "[stats,labeling]\n"
    output += "Par, Description,     Residues, Pairs\n"
    output += "N,   Total number,    {:>8}, {:>5}\n".format(stats["Na"], stats["Np"])
    output += "LI,  Invisible,       {:>8}, {:>5}\n".format(stats["LIa"], stats["LIp"])
    output += "LU,  Unique code,     {:>8}, {:>5}\n".format(stats["LUa"], stats["LUp"])
    output += "LN2, AA type of both, {:>8}, {:>5}\n".format(stats["LN2a"], stats["LN2p"])
    output += "LN1, AA type of last, {:>8}, {:>5}\n".format(stats["LN1a"], stats["LN1p"])
    output += "LA,  Ambiguous code,  {:>8}, {:>5}\n".format(stats["LAa"], stats["LAp"])
    return output


def get_stock_from_solution(solution):
    label_options = {}
    for residue in solution:
        labels = []
        pattern = solution[residue]
        for label_type in Constants.BASIC_TYPES:
            if label_type.name in pattern:
                labels.append(label_type)
        label_options[residue] = labels
    return label_options


def calc_price(prices_table, aa, pattern):
    price = 0
    for label in pattern:
        price += prices_table[aa][label]
    return price


def sort_residues(residue_list, residue_rank):
    # bubble sort for residues by residue rank.
    # used it because standard sorted() method gives randomized results
    residues = list(residue_list)
    for i in range(len(residues) - 1):
        for j in range(len(residues) - 1 - i):
            if residue_rank[i] < residue_rank[i + j + 1]:
                temp_res = residues[i]
                temp_rank = residue_rank[i]
                residue_rank[i] = residue_rank[i + j + 1]
                residues[i] = residues[i + j + 1]
                residue_rank[i + j + 1] = temp_rank
                residues[i + j + 1] = temp_res
    return residues, residue_rank
