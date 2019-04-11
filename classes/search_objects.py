import sys
import copy
import time
from .constants import Constants, Pattern, ELB
from scipy.optimize import linprog
from classes.ucsl_io import write_best_scheme, write_product_stats, write_products
from classes.sequence_specific import Residue2Label, calc_price

# add logging, especially in the reading of blocks
LOG_ITERATION = 10000


class Stock:

    def __init__(self, stock_table, prices_dict):
        self.stock_table = stock_table
        self.prices = prices_dict
        self.check_consistency()

    def check_consistency(self):
        if self.prices:
            for residue in self.stock_table:
                label_options = self.stock_table[residue]
                for label in label_options:
                    msg = "ERROR. Price is not specified for available '{}' label of '{}' residue type".format(label, residue)
                    try:
                        if self.prices[residue][label] < 0:
                            print(msg)
                            exit()
                    except KeyError:
                        print(msg)
                        exit()
        else:
            for residue in self.stock_table:
                label_options = self.stock_table[residue]
                curr_dict = {}
                for label in label_options:
                    curr_dict[label] = 0
                self.prices[residue] = curr_dict



class Scheme:

    def __init__(self, name, ncs, samples, patterns):
        self.name = name
        self.ncs = ncs
        self.codes = set()
        self.patterns = patterns
        self.samples = samples
        self.good = self.check_codes()
        self.simplified = {}
        self.simplified_str = ""
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
                if Pattern.pattern_bigger(self.patterns[i], self.patterns[i+j+1]):
                    temp_pattern = self.patterns[i]
                    self.patterns[i] = self.patterns[i+j+1]
                    self.patterns[i+j+1] = temp_pattern
        self.simplify()

    def is_subset_of(self, other_simple):
        for pattern in self.simplified:
            if pattern not in other_simple:
                return False
            if self.simplified[pattern] > other_simple[pattern]:
                return False
        return True

    def add_pattern(self, new_pattern):
        self.patterns.append(new_pattern)
        self.add_new_codes(new_pattern)
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

    def add_new_codes(self, new_pattern):
        for pattern in self.patterns:
            self.codes.add(self.ncs.calc_code(pattern, new_pattern))
            self.codes.add(self.ncs.calc_code(new_pattern, pattern))
        self.codes.add(self.ncs.calc_code(new_pattern, new_pattern))

    def add_pattern_list(self, pattern_list):
        for pattern in pattern_list:
            self.add_pattern(pattern)

    def simplify(self):
        self.simplified, self.simplified_str = \
            Pattern.simplify_list_of_patterns(self.patterns)

    def __hash__(self):
        return(hash(self.simplified_str))

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
        self.simplified_str = ""
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
        return "\n".join(self.patterns)

    def full_str(self):
        return "[ELB samples = {} patterns = {}]\n".format(self.samples, len(self.patterns)) \
                 + str(self)

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
        some_block = self.blocks[self.product_list[0][0]][self.product_list[0][1]][0]
        self.ncs_name = some_block.ncs_name
        self.deuterated = some_block.deuterated
        self.last_blocks = []
        self.product_samples = 0
        self.product_patterns = 1
        for i in range(len(self.product_list)):
            self.product_samples += self.product_list[i][0]
            self.product_patterns *= self.product_list[i][1]

    def __len__(self):
        length = 1
        for counter in self.max_counters:
            length *= counter
        return length

    def __next__(self):
        if self.overflow:
            raise StopIteration
        scheme = ELB([''], self.ncs_name, self.deuterated)
        self.last_blocks = []
        for i in range(len(self.counters)):
            j = len(self.counters) - i - 1
            samples_n = self.product_list[j][0]
            patterns_n = self.product_list[j][1]
            curr_block = self.blocks[samples_n][patterns_n][self.counters[j]]
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
         return ".".join("{}x{}".format(i, j) for i, j in self.product_list)


class BestScheme:

    def __init__(self, scheme, blocks, price, label_dict, residues):
        self.scheme = scheme
        self.blocks = blocks
        self.price = price
        self.label_dict = label_dict
        self.samples = scheme.samples
        self.residues = residues


class ProductFinder:

    def __init__(self, blocks, max_samples, min_patterns, equal=False):
        self.all_blocks = blocks
        self.max_samples = max_samples
        self.min_patterns = min_patterns
        self.equal = equal
        self.products = []
        self.product_lists = []
        self.block_types = []
        self.products_found = False
        self.make_block_types_list()

    def find_products(self):
        self.product_lists = []
        self.add_block_to_product([], 1, 0, self.max_samples)  # recursive search for all possible products
        for list_of_block_types in self.product_lists:
            self.products.append(Product(self.all_blocks, list_of_block_types))
        return self.products

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

            if samples == max_samples and ((self.equal and size == self.min_patterns)
                                           or (not self.equal and size >= self.min_patterns)):
                self.product_lists.append(new_list_of_types)
            else:
                self.add_block_to_product(new_list_of_types, size, samples, max_samples)

    def block_types_in_order(self, type_1, type_2):
        if type_2[0] > type_1[0]:
            return True
        elif type_2[0] == type_1[0]:
            if type_2[1] >= type_1[1]:
                return True
        return False

    def make_block_types_list(self):
        self.block_types = []
        for samples_n in self.all_blocks:
            pattern_num_list = []
            for patterns_n in self.all_blocks[samples_n]:
                pattern_num_list.append(patterns_n)
            pattern_num_list.sort(reverse=True)
            for patterns_n in pattern_num_list:
                self.block_types.append((samples_n, patterns_n))


class SchemeOptimizer:

    def __init__(self, task_params, logger):
        self.ncs = task_params["ncs"]
        self.residues = task_params["residues"]
        self.all_blocks = task_params["blocks"]
        self.price_table = task_params["prices"]
        self.job_name = task_params["jobname"]
        self.logger = logger
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
        self.first_products = True
        # self.cpu_count = multiprocessing.cpu_count()

    def scheme_checked(self, scheme, checked_schemes):
        for checked_scheme in checked_schemes:
            if scheme.simplified == checked_scheme \
                    or Pattern.first_scheme_subset(scheme.simplified, checked_scheme):
                return True
        return False

    def find_best_scheme(self):
        checked_schemes = []
        checked_number = 0
        schemes = 1
        scheme_found = False
        best_scheme = None
        price_optimizer = PriceOptimizer(self.ncs, self.price_table, self.residues)
        for product in self.products:
            product_name = "{}x{}= {}".format(product.product_samples, product.product_patterns, str(product))
            output = ("-"*70+"\nChecking product {}, total {} schemes\n"+"-"*70).format(
                product_name, len(product))
            self.logger.info(output)
            for scheme in product:
                curr_blocks = product.last_blocks
                output = "Checking scheme {:>5}/{:<5} - {:22} - ".format(schemes, self.schemes_total, product_name)
                status = ""

                if not self.scheme_checked(scheme, checked_schemes):
                    price_optimizer.minimize_price(scheme, curr_blocks)
                    status = "price calculated"
                    if not price_optimizer.success:
                        output += "STATUS: FAILED - Scheme can not be optimized"
                        status = "cannot be optimized"
                    elif not scheme_found or best_scheme.price > price_optimizer.best_scheme.price:
                        best_scheme = price_optimizer.best_scheme
                        scheme_found = True
                        output += "STATUS: {:12.5f} - New best price".format(best_scheme.price)
                    else:
                        output += "STATUS: {:12.5f}".format(price_optimizer.best_scheme.price)
                    checked_schemes.append(scheme.simplified)
                    checked_number += 1
                else:
                    output += "STATUS: SKIP - Equivalent scheme was already checked"
                    status = "equivalent"
                self.update_stats(scheme.samples, str(product.product_list), status)

                self.logger.info(output)
                output = "{}\nScheme consists of following blocks:\n".format(str(scheme))
                for block in curr_blocks:
                    output += str(block) + "\n"
                output += "---------"
                self.logger.debug(output)
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
        self.scheme_optimized = False
        max_samples = 1
        products_filename = "{}_products.txt".format(self.job_name)
        scheme_filename = "{}_scheme.txt".format(self.job_name)
        stats_filename = "{}_product_stats.txt".format(self.job_name)
        while not self.scheme_optimized:
            self.logger.info("Searching for products in max_samples={}".format(max_samples))
            self.find_products(max_samples)
            if self.products_found:
                mode = 'a'
                if self.first_products:
                    mode = 'w'
                    self.first_products = False
                output = write_products(self.products, max_samples, products_filename, mode=mode)
                self.logger.debug(output)
                self.logger.info("Schemes optimization in max_samples={}".format(max_samples))
                self.find_best_scheme()
            max_samples += 1
        output = write_best_scheme(self.best_scheme, scheme_filename)
        self.logger.debug(output)
        output = write_product_stats(self.stats, stats_filename)
        self.logger.debug(output)

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
                price = calc_price(self.prices, aa, pattern)
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
                if Pattern.simplify_pattern(pattern) == num_pattern:
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


class BlockFinder:

    def __init__(self, samples, ncs, min_depth, logger, elb_logger, block_finder_mode=False, min_t_free = -1):
        self.samples = samples
        self.block_finder_mode = block_finder_mode
        self.min_t_free = min_t_free
        self.check_t_free = False
        if min_t_free >=0:
            self.check_t_free = True
        self.index_of_type_T = Pattern.index_of_type(Constants.typeT)
        self.ncs = ncs
        self.min_depth = min_depth
        if self.min_depth <= 1:
            self.min_depth = 2
        self.scheme = Scheme("1", self.ncs, samples, [])
        self.patterns = [self.generate_patterns(self.samples)]
        self.depth = 0
        self.counter = [0]
        self.back_up_schemes = []
        self.result = {}
        self.results_found = 0
        self.output = ''
        self.iterator = 0
        self.max_depth = 0
        self.logger = logger
        self.elb_logger = elb_logger
        self.timer = time.time()

    def find(self):
        self.start_blockfinder()
        self.main_cycle()
        self.blockfinder_finished()

    def start_blockfinder(self):
        self.timer = time.time()

        if self.check_t_free:
            out = "[BlockFinder{}] started new search in {} samples with min_depth={}, min_t_free = {}".format(
                  self.samples, self.samples, self.min_depth, self.min_t_free)
        else:
            out = "[BlockFinder{}] started new search in {} samples with min_depth={}".format(
                  self.samples, self.samples, self.min_depth)
        self.logger.info(out)
        out = "[BlockFinder{}] total number of patterns is {}".format(
                  self.samples, len(self.patterns[0]))
        self.logger.info(out)
        self.logger.debug("List of all patterns used in all blocks:")
        for p in self.patterns[0]:
           self.logger.debug(p)

    def main_cycle(self):
        while True:
            self.next_iteration_output()
            patterns = self.patterns[self.depth]
            if self.depth == 0 and self.counter[0] + self.min_depth > len(patterns):
                break

            start_point = 1 + self.counter[self.depth]
            patterns_left = len(patterns) - start_point
            self.back_up_schemes.append(self.scheme.copy())
            self.scheme.add_pattern(patterns[self.counter[self.depth]])

            if patterns_left < self.min_depth - self.depth - 1:
                self.go_back()
                continue

            if patterns_left == 0:
                if len(self.scheme.patterns) >= self.min_depth:
                    self.save_result()
                self.go_back()
                continue

            next_patterns = self.get_next_patterns(patterns, patterns_left, start_point)

            flag_t_free = True
            if self.check_t_free:
                flag_t_free = self.check_have_enought_t_free(self.scheme, next_patterns)

            if next_patterns and flag_t_free:
                self.go_deeper(next_patterns)
            else:
                if len(self.scheme.patterns) >= self.min_depth:
                    self.save_result()
                self.go_parallel()


            self.check_max_depth()

    def next_iteration_output(self):
        self.iterator += 1
        if self.iterator % LOG_ITERATION == 0:
            out = "[BlockFinder{}] {:>9} {:>6d} sec ".format(self.samples, self.iterator,
                                                             int(time.time() - self.timer))
            out += "max_P={:<2} ELB_found= {:<6} ".format(self.max_depth + 1, self.results_found)
            for d in range(self.depth):
                out += " {:>3}/{:<3}".format(self.counter[d], len(self.patterns[d]) - self.min_depth + 1 + d)
            self.logger.info(out)

    def check_max_depth(self):
        if self.depth > self.max_depth:
            self.max_depth = self.depth
            if self.block_finder_mode:
                out = "[BlockFinder{}] New max depth: {}".format(self.samples, self.max_depth)
                self.logger.info(out)

    def check_have_enought_t_free(self, scheme, patterns_left):
        scheme_t, scheme_t_free = Pattern.count_type_in_list_of_simplified(scheme.simplified, self.index_of_type_T)
        left_t, left_t_free = Pattern.count_type_in_list_of_patterns(patterns_left, Constants.typeT)
        flag = scheme_t_free + left_t_free >= self.min_t_free
        return flag

    def get_next_patterns(self, patterns, patterns_left, start_point):
        next_patterns = []
        for i in range(patterns_left):
            if self.scheme.try_pattern(patterns[i + start_point]):
                next_patterns.append(patterns[i + start_point])
        return next_patterns

    def go_parallel(self):
        self.scheme = self.back_up_schemes[self.depth].copy()
        self.back_up_schemes.pop()
        self.counter[self.depth] += 1

    def go_back(self):
        self.depth -= 1
        self.patterns.pop()
        self.counter.pop()
        self.counter[-1] += 1
        self.back_up_schemes.pop()
        self.scheme = self.back_up_schemes[-1].copy()
        self.back_up_schemes.pop()

    def go_deeper(self, next_patterns):
        self.patterns.append(next_patterns)
        self.counter.append(0)
        self.depth += 1

    def blockfinder_finished(self):
        out = "[BlockFinder{}] finished search in {} samples after {:f} sec and {} iterations, {} ELB schemes found"
        out = out.format(self.samples, self.samples, time.time() - self.timer, self.iterator, self.results_found)
        self.logger.info(out)

    def generate_patterns(self, samples, top=True):
        if samples == 0:
            new_set = [""]
            return new_set
        current_set = self.generate_patterns(samples - 1, top=False)
        new_set = []
        for item in current_set:
            for option in self.ncs.label_types:
                new_pattern = item + option.name
                if top:
                    if self.ncs.check_power(new_pattern, self.min_depth):
                        new_set.append(new_pattern)
                else:
                    new_set.append(new_pattern)
        return new_set

    def save_result(self):
        if self.check_t_free and not self.check_have_enought_t_free(self.scheme, []):
            return
        depth_of_scheme = self.scheme.size()
        new_scheme = self.scheme.copy()
        new_scheme.sort()
        if depth_of_scheme in self.result:
            if new_scheme not in self.result[depth_of_scheme]:
                self.result[depth_of_scheme].update({new_scheme : 1 })
                self.write_result(new_scheme)
        else:
            self.result.update({depth_of_scheme: {new_scheme : 1}})
            self.write_result(new_scheme)

    def write_result(self, new_scheme):
        self.results_found += 1
        output = new_scheme.full_str()
        # self.logger.debug(output)
        self.elb_logger.info(output + "\n")

    # def write_result(self):
    #     out = "FindBlocks: finished after {} iterations\n".format(self.iterator)
    #     out += "FindBlocks: Evaluation time was {} seconds\n".format(time.time() - self.timer)
    #     out += "FindBlocks: Date/time is {}\n".format(time.strftime("%d-%m-%Y %H:%M:%S", time.gmtime()))
    #     self.logger.info(out)


class CLSequence:

    def __init__(self, seq):
        self.sequence = seq.upper()
        self.assignment = set()
        self.residues = dict()
        self.other_N = []
        self.other_C = []

    @property
    def residue_types_used(self):
        res_types = []
        for residue in Constants.RES_TYPES_LIST:
            if residue in self.sequence:
                res_types.append(residue)
        return res_types

    def sequence_stats(self, stock, assignment=set()):
        self.assignment = assignment
        self.residues = dict()
        self.other_N = []
        self.other_C = []

        for i in range(len(self.sequence) - 1):
            if (i + 1) in self.assignment:
                continue
            pair = self.sequence[i:i + 2]
            first_res = pair[0]
            second_res = pair[1]
            prices = {}
            if first_res in self.residues:
                self.residues[first_res].add_residue_after(second_res)
            else:
                label_options = stock.stock_table[first_res]

                if stock.prices:
                    prices = stock.prices[first_res]
                residue_obj = Residue2Label(first_res, label_options, prices)
                residue_obj.add_residue_after(second_res)
                self.residues[first_res] = residue_obj

            if second_res in self.residues:
                self.residues[second_res].add_residue_before(first_res)
            else:
                label_options = stock.stock_table[second_res]

                if stock.prices:
                    prices = stock.prices[second_res]
                residue_obj = Residue2Label(second_res, label_options, prices)

                residue_obj.add_residue_before(first_res)
                self.residues[second_res] = residue_obj

        # deal with Other N:
        for res in self.residues:
            res_obj = self.residues[res]
            if not res_obj.has_15n:
                for res2 in self.residues:
                    res_obj2 = self.residues[res2]
                    res_obj2.residues_after.discard(res)
                self.other_N.append(res)

        for res in self.residues:
            res_obj = self.residues[res]
            if not res_obj.has_13c:
                for res2 in self.residues:
                    res_obj2 = self.residues[res2]
                    res_obj2.residues_before.discard(res)
                    res_obj2.residues_before.add("Other")
                self.other_C.append(res)

        for res in self.residues:
            res_obj = self.residues[res]
            if res_obj.has_15n and res in res_obj.residues_before \
                    and "Other" in res_obj.residues_before and res_obj.has_double_label:
                res_obj.include_other = False

    def get_prices(self):
        prices_dict = {}
        for res_name in self.residues:
            res_obj = self.residues[res_name]
            prices_dict[res_name] = res_obj.labeling_prices
        return prices_dict


