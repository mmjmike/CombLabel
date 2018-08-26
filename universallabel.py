#!/usr/bin/python3

import copy
import sys
import datetime
import argparse
from scipy.optimize import linprog
from cl_errors import errors as err
from classes.constants import NCS, Spectrum, LabelType, Constants, Pattern
from classes.ucsl_io import Outputer, TaskReader
from classes.search_objects import Scheme, Product, BestScheme
import time
import concurrent.futures
import multiprocessing
import random


VERSION = "1.0.1"
CITATION = '''
#############################################################################
#                                                                           #
#        UCSL (Version {})- Calculation of Universal                     #
#        Combinatorial Labeling Schemes for Fast NMR Protein                #
#        Assignment.                                                        #
#                                                                           #
#        In case of using for scientific purpose, please cite our article:  #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#############################################################################
'''.format(VERSION)



class BlockFinder:

    def __init__(self, types, samples, ncs, min_depth, block_finder_mode, outputer=None):
        self.types = types
        # self.types_in_order()
        self.samples = samples
        self.ncs = ncs
        self.min_depth = min_depth
        self.scheme = Scheme("1", self.ncs, samples, [])
        self.patterns = [self.generate_patterns(self.samples)]
        # random.shuffle(self.patterns)
        # print(self.patterns)
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
        if self.min_depth == 1:
            self.min_depth = 2
        self.timer = time.time()


        out = "[BlockFinder{}] started new search in {} samples with min_depth={}\n".format(self.samples, self.samples, self.min_depth)
        self.outputer.write_data(out, files="lc")


        while True:
            self.iterator += 1

            if self.iterator % 10000 == 0:
                out = "[BlockFinder{}] {:>9} {:>6d} sec ".format(self.samples, self.iterator, int(time.time()-self.timer))
                out+= "max_P={:<2} ELB_found= {:<6} ".format(self.max_depth+1, self.results_found)
                for d in range(self.depth):
                    out += " {:>3}/{:<3}".format(self.counter[d], len(self.patterns[d])-self.min_depth+1+d)
                self.outputer.write_data(out, files="lc", timer=False)

            patterns = self.patterns[self.depth]
            if self.depth == 0 and self.counter[0] + self.min_depth > len(patterns):
                break
            self.back_up_schemes.append(self.scheme.copy())
            self.scheme.add_pattern(patterns[self.counter[self.depth]])
            if len(self.scheme.patterns) >= self.min_depth:
                self.save_result()
            next_patterns = []
            start_point = 1 + self.counter[self.depth]
            patterns_left = len(patterns) - start_point


            if patterns_left == 0 or patterns_left < self.min_depth - self.depth - 1:
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
                    out = "[BlockFinder{}] New max depth: {}".format(self.samples, self.max_depth)
                    self.outputer.write_data(out, files="lc")
        
        out = "[BlockFinder{}] finished search in {} samples after {} sec and {} iterations, {} ELB schemes found\n"
        out=out.format(self.samples, self.samples, int(time.time()-self.timer), self.iterator, self.results_found)
        self.outputer.write_data(out, files="lc")
        self.write_result()

    def generate_patterns(self, samples, top=True):
        if samples == 0:
            new_set = [""]
            return new_set
        current_set = self.generate_patterns(samples - 1, top=False)
        new_set = []
        for item in current_set:
            for option in self.types:
                new_pattern = item + option.name
                if top:
                    if self.ncs.check_power(new_pattern, self.min_depth):
                        new_set.append(new_pattern)
                else:
                    new_set.append(new_pattern)
        return new_set

    def save_result(self):
        depth_of_scheme = self.scheme.size()
        new_scheme = self.scheme.copy()
        new_scheme.sort()
        if depth_of_scheme in self.result:
            equal = False
            for scheme in self.result[depth_of_scheme]:
                if new_scheme == scheme:
                    equal = True
            if not equal:
                self.result[depth_of_scheme].append(new_scheme)
                self.results_found += 1
                output = "[ELB samples = {} patterns = {}]\n".format(new_scheme.samples, len(new_scheme.patterns))\
                         + str(new_scheme) + "\n"
                self.outputer.write_data(output, files="e")
        else:
            self.result.update({depth_of_scheme: [new_scheme]})
            self.results_found += 1
            output = "[ELB samples = {} patterns = {}]".format(new_scheme.samples, len(new_scheme.patterns)) \
                     + str(new_scheme) + "\n"
            self.outputer.write_data(output, files="e")

    def write_result(self):
        if self.block_finder_mode:
            out =  "FindBlocks: finished after {} iterations\n".format(self.iterator)
            out += "FindBlocks: Evaluation time was {} seconds\n".format(time.time()-self.timer)
            out += "FindBlocks: Date/time is {}\n".format(time.strftime("%d-%m-%Y %H:%M:%S", time.gmtime()))
            self.outputer.write_data(out, files="lc")
        # self.output = ''
        #
        # for depth in self.result:
        #     for scheme in self.result[depth]:
        #         # print("Output:")
        #         scheme.sort()
        #         self.output += '{} {}\n'.format(depth, scheme)
        #         if BLOCK_FIND:



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

    def uncode_scheme(self, x_coords, ):
        self.scheme.sort()
        primary_dict = {}
        label_dict = {}
        residues = []

        for i in range(len(x_coords)):
            if self.result.x[i] == 1:
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
        self.best_scheme = BestScheme(self.scheme, self.blocks, self.result.fun, label_dict, residues)

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
            label_type = Constants.TYPES[i]
            if aa == "P":
                label_type = Constants.PROLINE_SUBSTITUTE[label_type]
            curr_price = 0
            if number:
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
            self.three_letter_residues.append(Constants.TO_THREE_LETTER_CODE[res])
        self.three_letter_residues.sort()
        self.aa_number = len(self.residues)

        self.results = {}
        self.stats = {}
        self.block_types = []
        self.products = []
        self.all_schemes = []
        self.best_scheme = None
        self.schemes_total = 0
        self.scheme_optimized = False
        self.products_found = False
        self.first_output = True
        self.cpu_count = multiprocessing.cpu_count()

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
        schemes = 1
        scheme_found = False
        best_scheme = None
        price_optimizer = PriceOptimizer(self.ncs, self.price_table, self.residues)
        for product in self.products:
            output = "------------\nChecking schemes for product {}. Total {} schemes\n------------".format(str(product), len(product))
            self.outputer.write_data(output, files="cl")
            for scheme in product:
                curr_blocks = product.last_blocks
                output = "Checking scheme {}/{} - ".format(schemes, self.schemes_total)
                status = ""

                if not self.scheme_checked(scheme, checked_schemes):
                    price_optimizer.minimize_price(scheme, curr_blocks)
                    status = "price calculated"
                    if not price_optimizer.success:
                        output += "Scheme can not be optimized"
                        status = "cannot be optimized"
                    elif not scheme_found or best_scheme.price > price_optimizer.best_scheme.price:
                        best_scheme = price_optimizer.best_scheme
                        scheme_found = True
                        output += "New best price: {}".format(best_scheme.price)
                    else:
                        output += "Price: {}".format(price_optimizer.best_scheme.price)
                    checked_schemes.append(scheme.simplified)
                    checked_number += 1
                else:
                    output += "Equivalent scheme was already checked"
                    status = "equivalent"
                self.update_stats(scheme.samples, str(product.product_list), status)


                self.outputer.write_data(output, files="cl")
                output = "{}Scheme consists of following blocks:\n".format(str(scheme))
                for block in curr_blocks:
                    output += str(block) + "\n"
                output += "---------"
                self.outputer.write_data(output, files="l")
                schemes += 1
        self.scheme_optimized = scheme_found
        self.best_scheme = best_scheme

    def update_stats(self, samples, product_type, status):
        if samples not in self.stats:
            self.stats.update({samples: {}})
        if product_type not in self.stats[samples]:
            zero_dict = {
                "price calculated": 0,
                "equivalent": 0,
                "cannot be optimized": 0
            }
            self.stats[samples].update({product_type: zero_dict})
        self.stats[samples][product_type][status] += 1

    def run(self):
        if self.block_finder_mode:
            block_finder = BlockFinder(self.ncs.label_types, self.block_samples, self.ncs, self.block_min_depth, self.block_finder_mode, outputer=self.outputer)
            block_finder.find()
            self.output_blockfinder(block_finder)
        else:
            if self.calculate_blocks:
                output = "[NCS = {}]\n".format(self.ncs.name)
                self.outputer.write_data(output, files="lce")
                self.outputer.write_data("Searching for elementary blocks(ELB)", files="lc")
                self.find_blocks()
                # self.outputer.write_blocks(self.all_blocks, self.ncs)
            elif not self.only_blocks:
                self.outputer.write_block_stats(self.all_blocks)
                self.outputer.write_data("Clearing redundant blocks(ELB)", files="lc")
                self.clear_redundant_blocks()
                self.outputer.write_block_stats(self.all_blocks)
                self.outputer.write_data("Clearing product blocks(ELB)", files="lc")
                self.clear_product_blocks()
                self.outputer.write_block_stats(self.all_blocks)
                self.outputer.write_blocks(self.all_blocks, self.ncs)
            if self.only_blocks:
                return
            self.scheme_optimized = False
            max_samples = 1
            while not self.scheme_optimized:
                self.outputer.write_data("Searching for products in max_samples={}".format(max_samples), files="lc")
                self.find_products(max_samples)
                if self.products_found:
                    self.outputer.write_products(self.products, max_samples)
                    self.outputer.write_data("Schemes optimization in max_samples={}".format(max_samples), files="lc")
                    self.find_best_scheme()
                max_samples += 1
            self.outputer.write_best_scheme(self.best_scheme)
            self.outputer.write_product_stats(self.stats)
        self.finish_output()

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

    def find_products(self, max_samples):
        self.products_found = False
        prod_finder = ProductFinder(self.all_blocks, max_samples, self.aa_number)
        self.products = prod_finder.find_products()
        if self.products:
            self.products_found = True
        self.schemes_total = sum([len(product) for product in self.products])

    def find_blocks(self):
        # self.output("Search for elementary blocks with max {} samples".format(self.max_block_size))
        min_depth = 1
        self.all_blocks = {}
        for i in range(self.max_block_size):
            # self.output("Search in {} samples:".format(i+1))
            # self.output("Time in find_blocks: {} seconds\n".format(str(time.time()-t0)))
            block_finder = BlockFinder(self.ncs.label_types, i+1, self.ncs, min_depth, self.block_finder_mode, outputer=self.outputer)
            block_finder.find()
            result = block_finder.result
            self.all_blocks.update({i+1: result})
            if result:
                min_depth = max(result) + 1
            else:
                min_depth += 1

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
                        if block.is_subset_of(good_block.simplified):
                            block_good = False
                            break
                    if block_good:
                        good_blocks.append(block)
                        new_block_list.append(block)
                if not new_block_list:
                    self.all_blocks[samples_num].pop(pattern_num)
                else:
                    self.all_blocks[samples_num][pattern_num] = new_block_list
        self.clear_empty_block_types()

    def clear_product_blocks(self):
        blocks_samples = list(self.all_blocks.keys())
        blocks_samples.sort(reverse=True)
        for samples_num in blocks_samples:
            patterns_numbers = list(self.all_blocks[samples_num].keys())
            patterns_numbers.sort(reverse=True)
            for patterns_num in patterns_numbers:
                good_blocks = []
                for block in self.all_blocks[samples_num][patterns_num]:
                    block_good = True
                    prod_finder = ProductFinder(self.all_blocks, samples_num, patterns_num, equal=True)
                    products = prod_finder.find_products()
                    for product in products:
                        if len(product.product_list) == 1:
                            continue
                        if not block_good:
                            break
                        for product_block in product:
                            if product_block == block:
                                block_good = False
                                break
                    if block_good:
                        good_blocks.append(block)
                self.all_blocks[samples_num][patterns_num] = good_blocks
        self.clear_empty_block_types()

    def clear_blocks(self):
        pass

    def clear_empty_block_types(self):
        bad_samples_keys = []
        for samples_num in self.all_blocks:
            bad_patt_keys = []
            for patterns_num in self.all_blocks[samples_num]:
                if not self.all_blocks[samples_num][patterns_num]:
                    bad_patt_keys.append(patterns_num)
            for patterns_num in bad_patt_keys:
                self.all_blocks[samples_num].pop(patterns_num)
            if self.all_blocks[samples_num] == {}:
                bad_samples_keys.append(samples_num)
        for samples_num in bad_samples_keys:
            self.all_blocks.pop(samples_num)

    def finish_output(self):
        output = datetime.datetime.now().strftime('-----\nCalculation finished at %Y-%m-%d %H:%M:%S!')
        self.outputer.write_data(output, files="cl")
        self.outputer.close_files()

    def output_blockfinder(self, blockfinder):
        pass


def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", "-c", help='Specify the config file', type=str)
    # parser.add_argument('--samples', '-s', help='Maximum number of samples to start with', type=int, choise=range(1,10))
    parser.add_argument('--ncs', '-n', help='NMR Coding System (NCS)', type=str)
    parser.add_argument('--max_block_size', '-m', help='Maximum size of elementary blocks (in samples)', type=int, choices=range(1,9))
    parser.add_argument('--job_name', '-j', help='Job name', type=str)
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



def main():

    args = read_args()
    task_reader = TaskReader(args)
    outputer = Outputer(task_reader.output_parameters)
    outputer.open_files()
    outputer.write_data(CITATION, files="cle")
    task = Task(task_reader.task_parameters, outputer)
    task.run()

if __name__ == "__main__":
    main()
