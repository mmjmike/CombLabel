import sys
import copy
import time
from .constants import Constants, Pattern
from scipy.optimize import linprog
from classes.ucsl_io import write_best_scheme, write_product_stats, write_products




class Scheme:

    def __init__(self, name, ncs, samples, patterns):
        self.name = name
        self.ncs = ncs
        self.codes = set()
        self.patterns = patterns
        self.samples = samples
        self.good = self.check_codes()
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

    def is_subset_of(self, other_simple):
        for pattern in self.simplified:
            if pattern not in other_simple:
                return False
            if self.simplified[pattern] > other_simple[pattern]:
                return False
        return True

    def add_pattern(self, new_pattern):
        if self.try_pattern(new_pattern):
            self.patterns.append(new_pattern)
            self.codes.update(self.new_codes)

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

    @property
    def simplified(self):
        simplified = {}
        for pattern in self.patterns:
            simple_pattern = Pattern.simplify_pattern(pattern)
            if simplified != {} and simple_pattern in simplified:
                simplified[simple_pattern] += 1
            else:
                simplified.update({simple_pattern: 1})
        return simplified

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


class ELB:

    def __init__(self, patterns, ncs_name, deuterated=False):
        self.patterns = patterns
        self.ncs_name = ncs_name
        self.deuterated = deuterated

    def __str__(self):
        return "\n".join(self.patterns)

    @property
    def samples(self):
        if self.patterns:
            return len(self.patterns[0])
        else:
            return 0

    def full_str(self):
        return "[ELB samples = {} patterns = {}]\n".format(self.samples, len(self.patterns)) \
                 + str(self)

    def __mul__(self, other):
        new_patterns = []
        for pattern_1 in self.patterns:
            for pattern_2 in other.patterns:
                new_patterns.append(pattern_1 + pattern_2)
        return ELB(new_patterns, self.ncs_name, self.deuterated)

    def __eq__(self, scheme):
        return self.simplified == scheme.simplified

    @property
    def simplified(self):
        simplified = {}
        for pattern in self.patterns:
            simple_pattern = simplify_pattern(pattern)
            if simplified != {} and simple_pattern in simplified:
                simplified[simple_pattern] += 1
            else:
                simplified.update({simple_pattern: 1})
        return simplified

    def sort(self):
        for i in range(len(self.patterns)-1):
            for j in range(len(self.patterns)-1-i):
                if pattern_bigger(self.patterns[i], self.patterns[i+j+1]):
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
            for patterns_n in self.all_blocks[samples_n]:
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
            output = "------------\nChecking schemes for product {}. Total {} schemes\n------------".format(str(product), len(product))
            self.logger.info(output)
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
                if simplify_pattern(pattern) == num_pattern:
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


class BlockFinder:

    def __init__(self, samples, ncs, min_depth, logger, elb_logger, block_finder_mode=False):
        self.samples = samples
        self.block_finder_mode = block_finder_mode
        self.ncs = ncs
        self.min_depth = min_depth
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
        if self.min_depth == 1:
            self.min_depth = 2
        self.timer = time.time()

        out = "[BlockFinder{}] started new search in {} samples with min_depth={}".format(
                  self.samples, self.samples, self.min_depth)
        self.logger.info(out)
        out = "[BlockFinder{}] total number of patterns is {}".format(
                  self.samples, len(self.patterns[0]))
        self.logger.info(out)
        for p in self.patterns[0]:
           self.logger.info(p)
        
        while True:
            self.iterator += 1

            if self.iterator % 10000 == 0:
                out = "[BlockFinder{}] {:>9} {:>6d} sec ".format(self.samples, self.iterator,
                                                                 int(time.time() - self.timer))
                out += "max_P={:<2} ELB_found= {:<6} ".format(self.max_depth + 1, self.results_found)
                for d in range(self.depth):
                    out += " {:>3}/{:<3}".format(self.counter[d], len(self.patterns[d]) - self.min_depth + 1 + d)
                self.logger.info(out)

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
                self.depth -= 1
                self.patterns.pop()
                self.counter.pop()
                self.counter[-1] += 1
                self.back_up_schemes.pop()
                self.scheme = self.back_up_schemes[-1].copy()
                self.back_up_schemes.pop()
                continue

            for i in range(patterns_left):
                if self.scheme.try_pattern(patterns[i + start_point]):
                    next_patterns.append(patterns[i + start_point])

            if next_patterns == []:
                self.scheme = self.back_up_schemes[self.depth].copy()
                self.back_up_schemes.pop()
                self.counter[self.depth] += 1
            else:
                self.patterns.append(next_patterns)
                self.counter.append(0)
                self.depth += 1

            if self.depth > self.max_depth:
                self.max_depth = self.depth
                if self.block_finder_mode:
                    out = "[BlockFinder{}] New max depth: {}".format(self.samples, self.max_depth)
                    self.logger.info(out)

        out = "[BlockFinder{}] finished search in {} samples after {} sec and {} iterations, {} ELB schemes found"
        out = out.format(self.samples, self.samples, int(time.time() - self.timer), self.iterator, self.results_found)
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
                self.write_result(new_scheme)
        else:
            self.result.update({depth_of_scheme: [new_scheme]})
            self.write_result(new_scheme)

    def write_result(self, new_scheme):
        self.results_found += 1
        output = new_scheme.full_str()
        self.logger.debug(output)
        self.elb_logger.info(output + "\n")

    # def write_result(self):
    #     out = "FindBlocks: finished after {} iterations\n".format(self.iterator)
    #     out += "FindBlocks: Evaluation time was {} seconds\n".format(time.time() - self.timer)
    #     out += "FindBlocks: Date/time is {}\n".format(time.strftime("%d-%m-%Y %H:%M:%S", time.gmtime()))
    #     self.logger.info(out)



