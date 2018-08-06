import argparse
from classes.search_objects import Scheme
from classes.constants import NCS


class BlockFinder:
    def __init__(self, samples, ncs, min_depth, block_finder_mode, outputer=None):
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
        self.timer = time.time()
        self.max_depth = 0
        self.block_finder_mode = block_finder_mode
        self.outputer = outputer

        # print(self.patterns)

    def find(self):
        if self.min_depth == 1:
            self.min_depth = 2
        self.timer = time.time()

        out = "[BlockFinder{}] started new search in {} samples with min_depth={}\n".format(self.samples, self.samples,
                                                                                            self.min_depth)
        self.outputer.write_data(out, files="lc")

        while True:
            self.iterator += 1

            if self.iterator % 10000 == 0:
                out = "[BlockFinder{}] {:>9} {:>6d} sec ".format(self.samples, self.iterator,
                                                                 int(time.time() - self.timer))
                out += "max_P={:<2} ELB_found= {:<6} ".format(self.max_depth + 1, self.results_found)
                for d in range(self.depth):
                    out += " {:>3}/{:<3}".format(self.counter[d], len(self.patterns[d]) - self.min_depth + 1 + d)
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
                if self.scheme.try_pattern(patterns[i + start_point]):
                    next_patterns.append(patterns[i + start_point])

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
        out = out.format(self.samples, self.samples, int(time.time() - self.timer), self.iterator, self.results_found)
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
                output = "[ELB samples = {} patterns = {}]\n".format(new_scheme.samples, len(new_scheme.patterns)) \
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
            out = "FindBlocks: finished after {} iterations\n".format(self.iterator)
            out += "FindBlocks: Evaluation time was {} seconds\n".format(time.time() - self.timer)
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


def find_ncs(ncs_name):
    pass


def read_ncs_file(filename):
    pass


def good_filename(filename):
    return len(filename)


def run_blockfinder_once():
    pass


def find_blocks(ncs, max_block_size):

    # self.output("Search for elementary blocks with max {} samples".format(self.max_block_size))
    min_depth = 1
    all_blocks = {}
    for i in range(max_block_size):
        # self.output("Search in {} samples:".format(i+1))
        # self.output("Time in find_blocks: {} seconds\n".format(str(time.time()-t0)))
        block_finder = BlockFinder(i+1, ncs, min_depth)
        block_finder.find()
        result = block_finder.result
        all_blocks.update({i+1: result})
        if result:
            min_depth = max(result) + 1
        else:
            min_depth += 1


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("ncs", "-n", help='Specify NCS (NMR Coding System)', type=str)
    parser.add_argument("samples", "-s",
                        help='Specify maximal number of samples (range 1-9) or exact number of samples (if run in scheme finding mode)',
                        type=int)
    parser.add_argument("--minpatterns", "-p",
                        help='Specify minimal number of patterns to run in scheme finding mode',
                        type=int)
    return parser.parse_args()


def check_args(args):
    samples = args.samples
    exact_patterns = False

    if samples < 1 or samples > 9:
        print("Number of samples '{}' is not in range 1-9".format(samples))
        return False
    min_patterns = 0

    try:
        min_patterns = args.minpatterns
        exact_patterns = True
    except AttributeError:
        pass

    ncs_found, ncs_file = find_ncs(args.ncs)
    if not ncs_found:
        print("NCS file not found")
        return False


def main():
    args = get_args()
    if not check_args(args):
        return

    ncs_good, ncs = read_ncs_file(ncs_file)
    while True:
        input_filename = input("Please input filename for BlockFinder results: ")
        if not good_filename(input_filename):
            print("'{}' can't be used as a filename. Please enter a valid one (0 < length < 256; can't use following symbols: )")
            if not ask_to_continue_input():
                return
    if exact_patterns:
        run_blockfinder_once()
    else:
        find_blocks()
    print("Finished successfully")


if __name__ == "__main__":
    main()