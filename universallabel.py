#!/usr/bin/python3

import copy
import sys
import argparse
from scipy.optimize import linprog
from cl_errors import errors as err
from classes.constants import NCS, Spectrum, LabelType, Constants
from classes.ucsl_io import Outputer, TaskReader
from classes.search_objects import Scheme, Product, BestScheme
import time
import concurrent.futures
import multiprocessing


PRICES_FILE = "BEST_PRICE_1H.csv"
FULL_OUTPUT = True

JOB_NAME = "NC2-8-8"

OUTPUT_FILE = JOB_NAME + ".txt"
WRITE_TO_FILE = True
WRITE_TO_CONSOLE = True
CONFIG_FILE = "UCSL.task"
RESULT_FILE = JOB_NAME + "_results.txt"
BLOCK_FIND = False





class BlockFinder:

    def __init__(self, types, samples, ncs, min_depth, block_finder_mode, outputer=None):
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
        self.iterator = 0
        self.timer=time.time()
        self.max_depth = 0
        self.block_finder_mode = block_finder_mode
        self.outputer = outputer

        # print(self.patterns)

    def find(self):
        self.timer = time.time()
        out = "BlockFinder.find() started at {}\n".format(time.strftime("%d-%m-%Y %H:%M:%S", time.gmtime()))
        out += "BlockFinder: samples={} min_depth={}\n".format(self.samples, self.min_depth)
        self.outputer.write_data(out, files="l")
        while True:
            self.iterator += 1
            if self.block_finder_mode and self.iterator % 10000 == 0:
                out = "{:>8} {:>6d} sec  depth={:<2}".format(self.iterator, int(time.time()-self.timer), self.depth)
                for d in range(self.depth):
                    out += " {:>3}/{:<3}".format(self.counter[d], len(self.patterns[d])-self.min_depth+1+d)
                    self.outputer.write_data(out, files="l")
            patterns = self.patterns[self.depth]
            if self.depth == 0 and self.counter[0] + self.min_depth + 1 > len(patterns):
                break
            self.back_up_schemes.append(self.scheme.copy())
            self.scheme.add_pattern(patterns[self.counter[self.depth]])
            next_patterns = []
            start_point = 1 + self.counter[self.depth]
            patterns_left = len(patterns) - start_point
            # print(patterns_left, self.depth)
            if patterns_left == 0 or patterns_left < self.min_depth - self.depth - 1:
                if len(self.scheme.patterns) >= self.min_depth:
                    self.save_result()
                    # self.scheme.output()
                self.depth -= 1
                self.patterns.pop()
                self.counter.pop()
                # if len(self.counter) == 0:
                #     print("Fail")
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
            if self.depth > self.max_depth:
                self.max_depth = self.depth
                if self.block_finder_mode:
                    out += "New max depth: {}".format(self.max_depth)
                    self.outputer.write_data(out, files="l")

        if self.block_finder_mode:
            out = "FindBlocks finished after {} iterations\n".format(self.iterator)
            out += "FindBlocks: Evaluation time was {} seconds\n".format(time.time()-self.timer)
            out += "FindBlocks: Date/time is {}\n".format(time.strftime("%d-%m-%Y %H:%M:%S", time.gmtime()))
            self.outputer.write_data(out, files="l")

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
                if self.scheme == scheme:
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
                if BLOCK_FIND:
                    print(self.output)


class PriceOptimizer:

    def __init__(self, ncs, prices_table, aa_list):
        self.prices = prices_table
        self.aa_list = aa_list
        self.ncs = ncs
        self.success = False
        self.scheme = None
        self.best_scheme = None
        self.result = None
        self.patterns = []
        self.blocks = []

    def minimize_price(self, scheme, blocks):
        self.blocks = blocks
        self.success = False
        self.scheme = scheme
        self.patterns = []
        b_eq = [1 for _ in self.aa_list]
        b_ub = []
        A_eq = [[] for _ in range(len(self.aa_list))]
        price_sheet = []
        c = []
        x_coords = []
        for key in scheme.simplified:
            self.patterns.append(key)
            b_ub.append(scheme.simplified[key])
        A_ub = [[] for i in range(len(self.patterns))]
        for i in range(len(self.aa_list)):
            aa = self.aa_list[i]
            aa_prices = []
            for j in range(len(self.patterns)):
                pattern = self.patterns[j]
                price = self.calc_price(aa, pattern)
                aa_prices.append(price)
                if price > 0:
                    for k in range(len(self.aa_list)):
                        A_eq[k].append(0)
                    for k in range(len(self.patterns)):
                        A_ub[k].append(0)
                    A_eq[i][-1] = 1
                    A_ub[j][-1] = 1
                    c.append(price)
                    x_coords.append((i, j))
            price_sheet.append(aa_prices)

        self.result = linprog(c, A_ub, b_ub, A_eq, b_eq)
        self.success = self.result.success
        if self.success:
            self.uncode_scheme(x_coords)

    def uncode_scheme(self, x_coords):
        self.scheme.sort()
        primary_dict = {}
        label_dict = {}
        residues = []

        for i in range(len(x)):
            if x[i] == 1:
                res = self.aa_list[x_coords[i][0]]
                res = Constants.TO_THREE_LETTER_CODE[res]
                residues.append(res)
                num_pattern = self.patterns[x_coords[i][1]]
                primary_dict.update({res: num_pattern})

        residues.sort()
        patterns = list(self.scheme.patterns)
        for res in residues:
            for i in range(len(patterns)):
                num_pattern = primary_dict[res]
                pattern = patterns[i]
                if self.scheme.simplify_pattern(pattern) == num_pattern:
                    if res == "Pro":
                        pattern = self.substitute_pro(pattern)
                    label_dict.update({res: pattern})
                    patterns.pop(i)
                    break
        self.best_scheme = BestScheme(self.scheme, self.blocks, self.result.fun, label_dict)

    def substitute_pro(self, pattern):
        new_pattern = ""
        for char in pattern:
            new_pattern += Constants.PROLINE_SUBSTITUTE[char]
        return new_pattern

    def calc_price(self, aa, pattern):
        price = 0
        for i in range(len(pattern)):
            label = pattern[i]
            number = int(label)
            label_type = self.ncs.TYPES[i]
            if aa == "P":
                label_type = Constants.PROLINE_SUBSTITUTE[label_type]
            curr_price = self.prices[aa][label_type] * number
            if curr_price < 0:
                return -1
            price += curr_price
        return price


class Task:

    def __init__(self, task_params, outputer):
        self.ncs = task_params["NCS"]
        self.residues = task_params["aa_list"]
        self.max_block_size = task_params["max_block_size"]
        self.all_blocks = task_params["blocks"]
        self.only_blocks = task_params["only_blocks"]
        self.block_finder_mode = task_params["block_finder_mode"]
        self.calculate_blocks = task_params["calculate_blocks"]
        self.price_table = task_params["prices"]
        self.job_name = task_params["job_name"]
        self.block_samples = task_params["block_samples"]
        self.block_min_depth = task_params["block_min_depth"]

        self.outputer = outputer
        self.three_letter_residues = []
        for res in self.residues:
            self.three_letter_residues.append(self.ncs.to_three_letter_code[res])
        self.three_letter_residues.sort()
        self.aa_number = len(self.residues)

        self.results = {}
        self.block_types = []
        self.products = []
        self.product_lists = []
        self.all_schemes = []
        self.best_scheme = None

        self.product_schemes = 0
        self.scheme_optimized = False
        self.products_found = False
        self.first_output = True

        self.output_stream = open(OUTPUT_FILE, "w")
        self.logfile = open(JOB_NAME+".log", "w")
        self.cpu_count = multiprocessing.cpu_count()

    def __delete__(self):
        self.output_stream.close()

    def scheme_checked(self, scheme, checked_schemes):
        for checked_scheme in checked_schemes:
            if scheme.simplified == checked_scheme \
                    or self.first_scheme_subset(scheme.simplified, checked_scheme):
                return True
        return False

    def first_scheme_subset(self, scheme_1, scheme_2):
        for pattern in scheme_1:
            if pattern not in scheme_2:
                return False
            if scheme_1[pattern] > scheme_2[pattern]:
                return False
        return True

    def find_best_scheme(self):
        checked_schemes = []
        checked_number = 0
        scheme_found = False
        best_scheme = None
        price_optimizer = PriceOptimizer(self.ncs, self.price_table, self.residues)
        for product in self.products:
            for scheme in product:
                curr_blocks = product.last_blocks
                if not self.scheme_checked(scheme, checked_schemes):
                    price_optimizer.minimize_price(scheme, curr_blocks)
                    if not price_optimizer.success:
                        # signal that scheme can not be optimized
                        pass
                    elif not scheme_found or best_scheme.price > price_optimizer.best_scheme.price:
                        best_scheme = price_optimizer.best_scheme
                        scheme_found = True
                    checked_schemes.append(scheme.simplified)
                    checked_number += 1
        self.scheme_optimized = scheme_found
        self.best_scheme = best_scheme
        #     self.all_schemes = []
        #     self.obtain_scheme(product, [Scheme("", self.ncs, 0, [""])], [[]])
        #     self.output("\nChecking product: {}".format(product))
        #     self.output("total schemes for this product: {}".format(len(self.all_schemes)))
        #     for i in range(len(self.all_schemes)):
        #         scheme = self.all_schemes[i]
        #         scheme.sort()
        #         blocks = self.all_blocks[i]
        #         self.output("------------------------------------\n{} scheme. Consists of following blocks: ".format(checked_number+1))
        #         for block in blocks:
        #             self.output(str(block))
        #         self.output("\nScheme itself:")
        #         self.output(str(scheme))
        #         if scheme.simplified not in checked_schemes:
        #             optimal_scheme = price_optimizer.minimize_price(scheme)
        #             # if impossible to use due to lack in the stock, skip
        #             if not optimal_scheme[0].success:
        #                 self.output("Scheme can not be optimized")
        #             else:
        #                 label_dict = self.get_scheme_from_linprog(optimal_scheme, scheme)
        #                 self.output("\nScheme was successfully optimized:")
        #                 self.output_best_scheme(label_dict, scheme.samples, final=False)
        #                 self.output("Scheme price: {}".format(optimal_scheme[0].fun))
        #                 if not scheme_found:
        #                     best_scheme = optimal_scheme
        #                     best_scheme_patterns = copy.copy(scheme)
        #                     best_blocks = blocks
        #                     scheme_found = True
        #                 elif optimal_scheme[0].fun < best_scheme[0].fun:
        #                     best_scheme = optimal_scheme
        #                     best_scheme_patterns = copy.copy(scheme)
        #                     best_blocks = blocks
        #             checked_schemes.add(scheme.simplified)
        #         else:
        #             self.output("Equivalent scheme was already checked")
        #         checked_number += 1
        #
        #         self.output("{}/{} schemes checked".format(checked_number,
        #                                                    self.product_schemes))
        #         if scheme_found:
        #             self.output("Best price: {}".format(best_scheme[0].fun))
        #
        # if scheme_found:
        #     best_scheme_patterns.sort()
        #     label_dict = self.get_scheme_from_linprog(best_scheme, best_scheme_patterns)
        #     samples = best_scheme_patterns.samples
        #     self.output_best_scheme(label_dict, samples, best_blocks)
        # else:
        #     self.output("No schemes were found...")
        # self.output("________________________\nCalculation finished!")

    def run(self):
        if self.block_finder_mode:
            block_finder = BlockFinder(self.ncs.label_types, self.block_samples, self.ncs, self.block_min_depth)
            block_finder.find()
            self.output_blockfinder(block_finder)
        else:
            if self.calculate_blocks:
                self.find_blocks()
                self.clear_redundant_blocks()
                self.outputer.write_blocks(self.all_blocks)
            elif not self.only_blocks:
                self.read_blocks()
            if self.only_blocks:
                return
            self.output("Calculating all possible product types")

            self.scheme_optimized = False
            max_samples = 1
            while not self.scheme_optimized:
                self.find_products(max_samples)
                if self.products_found:
                    self.output_product_stats()
                    self.find_best_scheme()
                max_samples += 1
            self.outputer.write_best_scheme(self.best_scheme)

    # def get_scheme_from_linprog(self, best_scheme, scheme_patterns):
    #     simple_patterns = best_scheme[1]
    #     simple_patterns_num = len(simple_patterns)
    #     result = best_scheme[0]
    #     price_sheet = best_scheme[2]
    #     x = result.x
    #     scheme = {}
    #     patterns = list(scheme_patterns.patterns)
    #     skip = 0
    #     for i in range(self.aa_number):
    #         res = self.ncs.to_one_letter_code[self.three_letter_residues[i]]
    #         for j in range(simple_patterns_num):
    #             if price_sheet[i][j] < 0:
    #                 skip += 1
    #                 continue
    #             if x[i*simple_patterns_num+j-skip]:
    #                 simple_pattern = simple_patterns[j]
    #                 for k in range(len(patterns)):
    #                     pattern = patterns[k]
    #
    #                     if scheme_patterns.simplify_pattern(pattern) == simple_pattern:
    #                         if res == "P":
    #                             pattern = self.substitute_pro(pattern)
    #                         scheme[res] = ", ".join(list(pattern))
    #                         patterns.pop(k)
    #                         break
    #     return scheme

    def output_best_scheme(self, scheme, samples, blocks="", final=True):
        output = ""
        if final:
            output += '________________________\nBest scheme:\n[solution]'
        output += "\nRes"
        for i in range(samples):
            output += ",S{}".format(i + 1)
        output += "\n"
        for res in self.three_letter_residues:
            output += "{}, {}\n".format(res, scheme[self.ncs.to_one_letter_code[res]])
        if blocks:
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

    def add_block_to_product(self, curr_list_of_types, curr_size, curr_samples, max_samples):
        for block_type in self.block_types:
            if curr_list_of_types and not self.block_types_in_order(curr_list_of_types[-1], block_type):
                continue
            size = curr_size * block_type[1]
            samples = curr_samples + block_type[0]
            if samples > max_samples:
                continue
            new_list_of_types = copy.copy(curr_list_of_types)
            new_list_of_types.append(block_type)
            if size >= self.aa_number and samples == max_samples:
                self.product_lists.append(new_list_of_types)
            else:
                self.add_block_to_product(new_list_of_types, size, samples)

    def block_types_in_order(self, type_1, type_2):
        if type_2[0] > type_1[0]:
            return True
        elif type_2[0] == type_1[0]:
            if type_2[1] >= type_1[1]:
                return True
        return False

    def find_products(self, max_samples):
        self.products_found = False
        self.make_block_types_list()
        self.product_lists = []
        self.products = []
        self.add_block_to_product([], 1, 0, max_samples)  # recursive search for all possible products
        for list_of_block_types in self.product_lists:
            self.products.append(Product(self.all_blocks, list_of_block_types))
        if len(self.products) > 0:
            self.products_found = True

    def output_product_stats(self):
        self.product_schemes = 0
        for list_of_block_types in self.product_lists:
            number_of_schemes = 1
            for block_type in list_of_block_types:
                new_number = len(self.results[block_type[1]])
                number_of_schemes *= new_number
            self.product_schemes += number_of_schemes
            self.output(str(list_of_block_types) + ": {} scheme(s)".format(number_of_schemes))
        self.output("{} product types calculated:".format(len(self.products)))
        self.output("{} total labeling schemes to check".format(self.product_schemes))

    def make_block_types_list(self):
        self.block_types = []
        for samples_n in self.all_blocks:
            for patterns_n in self.all_blocks[samples_n]:
                self.block_types.append(tuple(samples_n, patterns_n))

    def find_blocks(self):
        # self.output("Search for elementary blocks with max {} samples".format(self.max_block_size))
        min_depth = 1
        self.all_blocks = {}
        for i in range(self.max_block_size):
            # self.output("Search in {} samples:".format(i+1))
            # self.output("Time in find_blocks: {} seconds\n".format(str(time.time()-t0)))
            block_finder = BlockFinder(self.ncs.label_types, i+1, self.ncs, min_depth)
            block_finder.find()
            result = block_finder.result
            self.all_blocks.update({i: result})
            if result:
                min_depth = max(result) + 1
            else:
                min_depth += 1


            # for key in result:
            #     self.types.append((i+1, key))
            #     blocks_total += len(result[key])

            # print(result)
        # sizes = []
        # for key in self.results:
        #     sizes.append(key)
        # sizes.sort()
        # for size in sizes:
        #     self.output("{} blocks {}x{} found".format(len(self.results[size]), self.results[size][0].samples, size))
        # self.output("{} elementary blocks found\n".format(blocks_total))

    def print_blocks(self):
        output = ''
        for key in self.results:
            output += str(key) + "\n"
            output += self.results[key] + "\n"
        self.output(output)

    def clear_redundant_blocks(self):
        for samples_num in self.all_blocks:
            patterns_numbers = list(self.all_blocks[samples_num].keys())
            patterns_numbers.sort(reverse=True)
            good_blocks = []
            for pattern_num in patterns_numbers:
                new_block_list = []
                for block in self.all_blocks[samples_num][pattern_num]:
                    block_good = True
                    for good_block in good_blocks:
                        if block.is_subset_of(good_block):
                            block_good = False
                            break
                    if block_good:
                        good_blocks.append(block)
                        new_block_list.append(block)
                if new_block_list == []:
                    self.all_blocks[samples_num].pop(pattern_num)
                else:
                    self.all_blocks[samples_num][pattern_num] = new_block_list

    def read_blocks(self):
        pass

    def output(self, text, result=False):
        mode = 'a'
        if self.first_output:
            mode = 'w'
            self.first_output = False
        if WRITE_TO_FILE:
            self.write_to_file(OUTPUT_FILE, text + "\n", mode=mode)
        if result:
            self.write_to_file(RESULT_FILE, text + "\n", mode=mode)
        if WRITE_TO_CONSOLE:
            print(text)

    def write_to_file(self, filename, output, mode='w'):

        #f = open(filename, mode)
        self.output_stream.write(output)
        self.output_stream.flush()
        #f.close()

    def output_blockfinder(self, blockfinder):
        pass





# class Stock:
#     '''
# Stock object
#
# Contains the amino acid stock
# and also prices for each type of labeling
# '''
#
#     def __init__(self, task):
#         self.task = task
#         self.label_dict = {}
#         self.label_options = {}
#         self.label_num_dict = {}
#         self.price_dict = {}
#         self.usage = {}
#         self.label_types = []
#         self.price_label_types = []
#         for res in self.task.res_types:
#             self.usage.update({res: 1})
#
#     def read_stock(self, lines):
#         d = []
#         first_entry = True
#         row_len = 0
#         for line in lines:
#             if line[0] != "#" and line != "":
#                 s = [x.strip() for x in line.split(",")]
#                 if first_entry:
#                     row_len = len(s)
#                     first_entry = False
#                 elif len(s) != row_len:
#                     message = "ERROR in stock file!\nThe lenths of rows are not equal"
#                     raise err.ReadConfigError(msg=message)
#                 d.append(s)
#
#         try:
#             self.label_types = d[0][1:]
#             for label_type in self.label_types:
#                 if label_type not in self.task.label_types:
#                     raise err.ReadConfigError(msg="ERROR! Wrong labeling type in stock file.")
#             for i in range(len(d) - 1):
#                 res_name = d[i + 1][0]
#                 if res_name in self.task.res_types:
#                     res = res_name
#                 elif res_name.upper() in [r.upper() for r in self.task.res_types_3]:
#                     res = self.task.to_one_letter_code[res_name]
#                 else:
#                     raise err.ReadConfigError(msg="ERROR! Wrong residue name in stock file.")
#                 self.label_dict.update({res: ''.join(d[i + 1][1:])})
#             self._generate_label_options()
#         except IndexError:
#             raise err.ReadConfigError(msg="ERROR in stock file.\nPlease check the length of each row.\nUse comma as separator")
#         except err.ReadConfigError as e:
#             raise err.ReadConfigError(msg=e.msg)
#
#     def read_from_table(self, table):
#         #self.label_types = ['X', 'N', 'C', 'D']
#         for i in range(len(self.task.res_types)):
#             res = self.task.res_types[i]
#             table_part = [table[j][i] for j in range(len(table))]
#             self.label_dict.update({res: ''.join(["1" if cell else "0" for cell in table_part])})
#         self._generate_label_options()
#
#     def read_prices_from_table(self, table):
#         pass
#
#     def _generate_label_options(self):
#         for residue in self.task.res_types:
#             if residue in self.label_dict:
#                 stock_to_change = self.label_dict[residue]
#                 option = []
#                 ## Labeling types are in preferred order,
#                 ##but due to success with this software the order doesn't matter.
#                 for label in self.task.label_types:
#                     if label in self.label_types and label in self.task.coding_table.label_types_list:
#                         label_index = self.label_types.index(label)
#                         if stock_to_change[label_index] == '1':
#                             option.append(label)
#                 label_num_dict = {}
#                 for i in range(len(option)):
#                     label = option[i]
#                     label_num_dict[label] = i
#                 self.label_options.update({residue: option})
#                 self.label_num_dict.update({residue: label_num_dict})
#
#     def read_prices(self, lines):
#         d = []
#         for line in lines:
#             a = line[0]
#             if line[0] != "#" and line != "":
#                 s = [x.strip() for x in line.split(",")]
#                 d.append(s)
#         self.price_label_types = d[0][1:]
#         for label in self.label_types:
#             if label not in self.price_label_types:
#                 raise err.ReadConfigError(msg="ERROR! Prices are not specified for '" + label + "' label type")
#
#         # ADD check if file format is correct
#         try:
#             for i in range(len(d) - 1):
#                 curr_dict = {}
#                 for j in range(len(self.price_label_types)):
#                     try:
#                         price = float(d[i + 1][j + 1])
#                     except ValueError:
#                         raise err.ReadConfigError(msg="ERROR! Price must be set in digits\nPlease check price file (row {}; col {})".format(i+2, j+2))
#                     curr_dict.update({self.price_label_types[j]: price})
#                 res_name = d[i + 1][0]
#                 if res_name in self.task.res_types:
#                     residue_type = res_name
#                 elif res_name.upper() in [r.upper() for r in self.task.res_types_3]:
#                     residue_type = self.task.to_one_letter_code[res_name]
#                 else:
#                     raise err.ReadConfigError(msg="ERROR! Wrong residue name in price file.")
#                 for label_type in self.label_options[residue_type]:
#                     if curr_dict[label_type] < 0:
#                         raise err.ReadConfigError(msg="ERROR! Price is not specified or negative for '"
#                               + residue_type + "'-residue's '"
#                               + label_type + "' label type")
#                 self.price_dict.update({residue_type: curr_dict})
#         except IndexError:
#             raise err.ReadConfigError(msg="ERROR in price file.\nPlease check the length of each row.\nUse comma as separator")
#         for res in self.label_options:
#             for option in self.label_options[res]:
#                 try:
#                     self.price_dict[res][option]
#                 except KeyError:
#                     raise err.ReadConfigError(msg="ERROR in price file.\n"
#                                                   "Price is not specified for {} label "
#                                                   "option of {}".format(option, self.task.to_three_letter_code[res]))
#
#     def write_to_file(self, filename, mode='w'):
#
#         #ADD check if file exists
#
#         f = open(filename, mode)
#         f.write('Res\t' + '\t'.join(self.task.label_types) + '\n')
#         for residue in self.task.res_types:
#             line = residue
#             for label in self.task.label_types:
#                 line += '\t'
#                 if residue in self.price_dict:
#                     line += str(self.price_dict[residue][label])
#                 else:
#                     line += '0'
#             line += '\n'
#             f.write(line)
#         f.flush()
#         f.close()
#
#     def log_print(self):
#         i = 1
#         for key in self.label_dict:
#             self.task.logfile.write (str(i) + "." + key + ":" + self.label_dict[key])
#             i += 1
#         self.task.logfile.write("\n")
#         self.task.logfile.flush()



def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", "-c", help='Specify the config file', type=str)
    # parser.add_argument('--samples', '-s', help='Maximum number of samples to start with', type=int, choise=range(1,10))
    parser.add_argument('--ncs', '-n', default='NC2', help='NMR Coding System (NCS)', type=str)
    parser.add_argument('--max_block_size', '-m', help='Maximum size of elementary blocks (in samples)', type=int, choise=range(1,9))
    parser.add_argument('--job_name', '-j', default='Default_task', help='Job name', type=str)
    parser.add_argument('--stock', '-t', help='Input stock file', type=str)
    parser.add_argument('--blocks', '-b', help='Input blocks file', type=str)
    parser.add_argument('--only_blocks', '-o', help='Program runs only block calculation without further schemes optimization', action="store_true")
    parser.add_argument('--verbose', '-v',
                        help='Verbose output to console',
                        action="store_true")
    parser.add_argument('--silent', '-s',
                        help='No output to console (overrides verbosity)',
                        action="store_true")
    parser.add_argument("--block_finder_mode", "-f", nargs=2, type=int,
                        help="Enables running in only block finder mode. Required two integer parameters: samples (1<=s<=9) and minimal depth of search (1<=m<=20)")
    parser.add_argument('--price_table', '-p', help='Input prices list file for amino acids', type=str)
    return parser.parse_args()


constants = Constants()
task1 = Task(NC2, RES_TYPES_LIST, 9, 5)

task2 = Task(NCD2, RES_TYPES_LIST, 7, 4)
task3 = Task(NCD4, RES_TYPES_LIST, 6, 3)
task4 = Task(NCD6, RES_TYPES_LIST, 5, 3)
task5 = Task(NCDA8, RES_TYPES_LIST, 4, 3)
task6 = Task(TSF12, RES_TYPES_LIST, 4, 2)
block_find = BlockFinder([typeX, typeN, typeC], 1, NC2, 1)
block_find8_20 = BlockFinder([typeN, typeC], 8, NC2, 20)
test_task = Task(NC2, RES_TYPES_LIST, 9, 5)


def main():

    args = read_args()
    task_reader = TaskReader(args)
    outputer = Outputer(task_reader.output_parameters)
    outputer.open_files()
    task = Task(task_reader.task_parameters)
    task.run()

if __name__ == "__main__":
    main()
