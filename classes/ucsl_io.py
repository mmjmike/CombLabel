from cl_errors import errors as err
from classes.constants import Constants
from classes.search_objects import Scheme
import re
import logging
import datetime
import random


class Outputer:

    def __init__(self, output_parameters):
        self.job_name = output_parameters["job_name"]
        self.verbose = output_parameters["verbose"]
        self.files = output_parameters["files"]
        self.silent = output_parameters["silent"]

        self.log_filename = self.job_name + ".log"
        self.blocks_filename = self.job_name + "_elb.txt"
        self.scheme_filename = self.job_name + "_scheme.txt"
        self.products_filename = self.job_name + "_products.txt"
        self.open_files()
        self.opened_files = {}

    def set_verbose(self):
        self.verbose = True

    def silent(self):
        self.silent = True

    def open_files(self):
        types = list(self.files)
        if "l" in types:
            self.opened_files.update({"l": open(self.log_filename, "w")})
        if "e" in types:
            self.opened_files.update({"e": open(self.blocks_filename, "w")})
        if "p" in types:
            self.opened_files.update({"p": open(self.products_filename, "w")})
        if "s" in types:
            self.opened_files.update({"s": open(self.scheme_filename, "w")})

    def write_data(self, output, files="l"):
        types = list(files)
        if "l" in types and "l" in self.opened_files:
            self.opened_files["l"].write(output)
            self.opened_files["l"].flush()
        if "e" in types and "e" in self.opened_files:
            self.opened_files["e"].write(output)
            self.opened_files["e"].flush()
        if "p" in types and "p" in self.opened_files:
            self.opened_files["p"].write(output)
            self.opened_files["p"].flush()
        if "s" in types and "s" in self.opened_files:
            self.opened_files["s"].write(output)
            self.opened_files["s"].flush()
        if "c" in types and self.verbose:
            print(output)

    def close_files(self):
        for key in self.opened_files:
            self.opened_files[key].close()

    def write_best_sheme(self, scheme):
        pass

    def write_products(self, products):
        pass

    def write_elb(self, blocks):
        pass




class TaskReader:

    def __init__(self, args):
        self.args = args
        self.config_exists = False
        self.config = ""
        try:
            self.config = self.args.config
            f = open(self.args.config)
        except AttributeError:
            pass
        except FileNotFoundError:
            pass
        else:
            self.config_exists = True
            f.close()
        self.silent = False
        self.verbose = False
        self.only_blocks = False
        self.block_finder_mode = False
        self.calculate_blocks = True
        self.block_samples = 1
        self.block_min_depth = 1
        self.max_block_size = 1
        self.blocks = {}
        self.prices = {}
        self.stock = {}
        self.mode = "optimize_scheme"
        self.job_name = ""
        self.aa_list = Constants.RES_TYPES_LIST
        flag_list = [False for _ in range(4)]
        flag_list += ["", ""]
        self.parameters_read = {
            "job_name": list(flag_list) + [self.check_job_name],
            "NCS": list(flag_list) + [self.check_ncs],
            "spectra": list(flag_list) + [self.check_spec_vec],
            "label_types": list(flag_list) + [self.check_types_vec],
            "aa_list": list(flag_list) + [self.check_aa_list],
            "price_table": list(flag_list) + [self.check_prices_file],
            "stock": list(flag_list) + [self.check_stock_file],
            "max_block_size": list(flag_list) + [self.check_block_size],
            "blocks": list(flag_list) + [self.check_blocks_file]
        }
        self.ncs_specified = False
        self.parameters_good = False
        self.ncs = Constants.NC2

        self.create_task()

    def check_job_name(self, job_name):
        return len(job_name) > 0

    def check_ncs(self, ncs):
        return ncs.upper() in Constants.NCS_NAMES

    def check_spec_vec(self, spec_vec):
        if len(spec_vec) < 2:
            return False
        for number in spec_vec:
            try:
                num = int(number)
                if num < 0:
                    return False
            except ValueError:
                return False
        return True

    def check_types_vec(self, types_vec):
        if len(types_vec) != 8:
            return False
        for number in types_vec:
            try:
                num = int(number)
                if num < 0:
                    return False
            except ValueError:
                return False
        return True

    def check_prices_file(self, prices_file):
        try:
            f = open(prices_file, "r")
            f.close()
        except FileNotFoundError:
            return False
        return self.read_prices(prices_file)[0]

    def read_prices(self, prices_file):
        msg = ""
        prices = {}
        lines = self.read_lines(prices_file)
        if len(lines) < 3:
            msg = "Prices file '{}' is too short".format(prices_file)
            return [False, msg, prices]
        d = []
        line_length = 0
        first_line = True
        for line in lines:
            s = [x.strip() for x in line.split(",")]
            if first_line:
                line_length = len(s)
                first_line = False
            else:
                if len(s) != line_length:
                    msg = "Not equal length of lines in prices file '{}'".format(prices_file)
                    return [False, msg, prices]
            d.append(s)
        price_label_types = d[0][1:]
        for l_type in price_label_types:
            if l_type.upper() not in Constants.TYPES_NAMES:
                msg = "Incorrect label type {} in prices file '{}'".format(l_type, prices_file)
                return [False, msg, prices]
        if len(price_label_types) < 2:
            msg = "Too few labeling types specified in prices file '{}'".format(prices_file)
            return [False, msg, prices]

        for i in range(len(d) - 1):
            curr_dict = {}
            for j in range(len(price_label_types)):
                try:
                    price = float(d[i + 1][j + 1])
                except ValueError:
                    msg="ERROR! Price must be set in digits"
                    msg += "\nPlease check price file '{}' (row {}; col {})".format(prices_file,
                            i + 2, j + 2)
                    return [False, msg, prices]
                curr_dict.update({price_label_types[j]: price})
            residue_type = d[i + 1][0]
            if residue_type not in Constants.RES_TYPES_LIST and residue_type not in Constants.RES_TYPES_THREE:
                msg = "Wrong residue type '{}' in prices file '{}'".format(residue_type, prices_file)
                return [False, msg, prices]
            self.prices.update({residue_type: curr_dict})
        return [True, msg, prices]

    def check_stock_file(self, stock_file):
        try:
            f = open(stock_file, "r")
            f.close()
        except FileNotFoundError:
            return False
        return self.read_stock(stock_file)[0]

    def read_stock(self, stock_file):
        msg = ""
        label_dict = {}
        lines = self.read_lines(stock_file)
        if len(lines) < 3:
            msg = "Stock file '{}' is too short".format(stock_file)
            return [False, msg, label_dict]
        d = []
        line_length = 0
        first_line = True
        for line in lines:
            s = [x.strip() for x in line.split(",")]
            if first_line:
                line_length = len(s)
                first_line = False
            else:
                if len(s) != line_length:
                    msg = "Not equal length of lines in stock file '{}'".format(stock_file)
                    return [False, msg, label_dict]
            d.append(s)
        stock_label_types = d[0][1:]

        for l_type in stock_label_types:
            if l_type.upper() not in Constants.TYPES_NAMES:
                msg = "Incorrect label type {} in stock file '{}'".format(l_type, stock_file)
                return [False, msg, label_dict]
        if len(stock_label_types) < 2:
            msg = "Too few labeling types specified in stock file '{}'".format(stock_file)
            return [False, msg, label_dict]

        for i in range(len(d) - 1):
            curr_list = []
            for j in range(len(stock_label_types)):
                try:
                    presence = int(d[i + 1][j + 1])
                    if presence == 1:
                        curr_list.append(stock_label_types[j])
                    elif presence == 0:
                        pass
                    else:
                        raise ValueError
                except ValueError:
                    msg="ERROR! Presence/absence of labeling type must be set in 1/0"
                    msg += "\nPlease check stock file '{}' (row {}; col {})".format(stock_file,
                            i + 2, j + 2)
                    return [False, msg, label_dict]
            residue_type = d[i + 1][0]
            if residue_type not in Constants.RES_TYPES_LIST and residue_type not in Constants.RES_TYPES_THREE:
                msg = "Wrong residue type '{}' in stock file '{}'".format(residue_type, stock_file)
                return [False, msg, label_dict]
            label_dict.update({residue_type: curr_list})
        return [True, msg, label_dict]

    def check_block_size(self, block_size):
        try:
            size = int(block_size)
        except ValueError:
            return False
        return 0 < size < 10

    def check_blocks_file(self, block_file):
        try:
            f = open(block_file, "r")
            f.close()
        except FileNotFoundError:
            return False
        return self.read_blocks(block_file)[0]

    def read_block_file(self, block_file):
        lines = self.read_lines(block_file)
        i = 0
        msg = ''
        blocks = {}
        blocks_num = 0
        NCS_found = False
        ncs = Constants.NC2
        ncs_regular = re.compile('\\[NCS = \\w+\\]')
        while i < len(lines):
            if not NCS_found:
                ncs_match = ncs_regular.match(lines[i])
                if str(ncs_match) == 'None':
                    i += 1
                    continue
                else:
                    ncs_name = ncs_match.group()[ncs_match.start()+7:-1].upper()
                    if ncs_name not in Constants.NCS_NAMES:
                        msg = "Warning! Wrong NCS '{}' is specified in blocks file {}".format(ncs_name, block_file)
                        continue
                    else:
                        NCS_found = True
                        ncs = [ncs_curr for ncs_curr in Constants.LIST_OF_NCS if ncs_curr.name == ncs_name][0]
                        i += 1
                        continue
            if NCS_found and re.search(r'\[ELB samples = \d patterns = \d+\]', lines[i]):
                blocks_num += 1
                numbers = re.findall(r'\d+', lines[i])
                samples_num = int(numbers[0])
                patterns_num = int(numbers[1])
                patterns = []
                good_block = True
                for j in range(patterns_num):
                    i += 1
                    pattern = lines[i].strip()
                    if len(pattern) != samples_num:
                        msg = "Warning! The number of samples in block {} doesn't match " \
                              "\nwith specified value in blocks file {}.".format(blocks_num, block_file)
                        good_block = False
                        # return [False, msg, blocks]
                if good_block:
                    block = Scheme("", ncs, samples_num, patterns)
                    if samples_num not in blocks:
                        blocks.update({samples_num: {patterns_num: [block]}})
                    else:
                        if patterns_num not in blocks[samples_num]:
                            blocks[samples_num].update({patterns_num: [block]})
                        else:
                            blocks[samples_num][patterns_num].append(block)
            else:
                i += 1
        return [NCS_found, msg, blocks]

    def create_task(self):
        if self.config_exists:
            lines = self.read_lines(self.args.config)
            self.read_config(lines)
        self.read_args()
        self.check_parameters()
        self.parameters_good, msg = self.check_consistency()
        if self.parameters_good:
            self.read_values()
        else:
            self.print_errors()

    def read_ncs(self):
        if self.parameters_read["NCS"][1]:
            self.ncs = [ncs for ncs in Constants.LIST_OF_NCS if self.parameters_read["NCS"][4] == ncs.name][0]
        if self.parameters_read["NCS"][3]:
            self.ncs = [ncs for ncs in Constants.LIST_OF_NCS if self.parameters_read["NCS"][5] == ncs.name][0]

    def read_job_name(self):
        if self.parameters_read["job_name"][1]:
            self.job_name = self.parameters_read["job_name"][4]
        if self.parameters_read["job_name"][3]:
            self.job_name = self.parameters_read["job_name"][5]
        if self.job_name == "":
            self.job_name = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')
            self.job_name += "_{:02}".format(random.randrange(0, 100))

    def set_max_block_size(self):
        if self.parameters_read["max_block_size"][1]:
            self.max_block_size = int(self.parameters_read["max_block_size"][4])
        if self.parameters_read["max_block_size"][3]:
            self.max_block_size = int(self.parameters_read["max_block_size"][5])

    def read_aa_list(self):
        if self.parameters_read["aa_list"][1]:
            self.aa_list = self.parameters_read["aa_list"][4].split()
        if self.parameters_read["aa_list"][3]:
            self.aa_list = self.parameters_read["aa_list"][5].split()

    def read_values(self):
        self.read_ncs()
        self.read_job_name()
        if not self.block_finder_mode:
            if self.calculate_blocks:
                self.set_max_block_size()
            else:
                if self.parameters_read["blocks"][1]:
                    blocks_file = self.parameters_read["blocks"][4]
                if self.parameters_read["blocks"][3]:
                    blocks_file = self.parameters_read["blocks"][5]
                self.blocks = self.read_blocks(blocks_file)[2]
            if not self.only_blocks:
                if self.parameters_read["price_table"][1]:
                    prices_file = self.parameters_read["price_table"][4]
                if self.parameters_read["price_table"][3]:
                    prices_file = self.parameters_read["price_table"][5]
                self.prices = self.read_prices(prices_file)[2]

                self.read_aa_list()

                stock_read = False
                if self.parameters_read["stock"][1]:
                    stock_file = self.parameters_read["stock"][4]
                    stock_read = True
                if self.parameters_read["stock"][3]:
                    stock_file = self.parameters_read["stock"][5]
                    stock_read = True
                if stock_read:
                    self.stock = self.read_stock(stock_file)[2]

    def read_config(self, lines):
        parameters = {}
        for line in lines:
            split_line = line.split()
            if len(split_line) > 1:
                parameters[split_line[0]] = split_line[1:]

        for key in self.parameters_read:
            if key in parameters:
                self.parameters_read[key][0] = True
                self.parameters_read[key][4] = parameters[key]

        #
        #
        # try:
        #     job_name = parameters["job_name"]
        #     self.job_name = job_name[0]
        # except (KeyError, IndexError):
        #     print("Warning! Job name is not specified in the config file")
        #
        # try:
        #     aa_list = parameters["aa_list"]
        #     if not self.check_aa_list(aa_list):
        #         print("Warning! Amino-acid types are not specified correctly in the config file")
        #     else:
        #         self.aa_list = aa_list
        # except KeyError:
        #     print("Warning! Amino-acid types are not specified in the config file")
        #
        # try:
        #     ncs_name = parameters["NCS"]
        #
        #     for ncs in Constants.LIST_OF_NCS:
        #         if ncs_name == ncs.name:
        #             self.ncs_specified = True
        #             self.ncs = ncs
        #             break
        #     if not self.ncs_specified:
        #         print("Warning! NCS is not specified correctly in the config file\n NC2 system taken as default value")
        # except KeyError:
        #     print("Warning! NCS is not specified in the config file")
        #
        # if not self.ncs_specified:
        #     try:
        #         spec_vec = parameters["spectra_types"]
        #         if spec_vec != []:
        #             if int(spec_vec[0]) == 0:
        #
        #
        #     except KeyError:
        #         print("Warning! Spectra types are not specified in the config file")
        #
        #         raise err.ReadConfigError(msg="ERROR! Spectra types are not specified in the config file")
        #     try:
        #         self.read_spec_vec(spec_vec)
        #     except err.ReadConfigError as e:
        #         raise err.ReadConfigError(msg=e.msg)
        #
        #     try:
        #         lab_vec = parameters["labeling_types"]
        #     except KeyError:
        #         raise err.ReadConfigError(msg="ERROR! Labeling types are not specified in the config file")
        #     if lab_vec == []:
        #         raise err.ReadConfigError(msg="ERROR! Labeling types are not specified in the config file")
        #     try:
        #         self.read_label_vec(lab_vec)
        #     except err.ReadConfigError as e:
        #         raise err.ReadConfigError(msg=e.msg)
        #
        #
        #
        # try:
        #     max_block_size = parameters["max_block_size"]
        # except KeyError:
        #     raise err.ReadConfigError(msg="ERROR! Max block size is not specified in the config file")
        # if max_block_size == []:
        #     raise err.ReadConfigError(msg="ERROR! Max block size is not specified in the config file")
        # try:
        #     self.set_block_size(max_block_size[0])
        # except err.ReadConfigError as e:
        #     raise err.ReadConfigError(msg=e.msg)
        #
        #
        #
        # try:
        #     self.make_coding_table()
        # except err.LabelPowerError as e:
        #     raise err.ReadConfigError(msg=e.msg)
        #
        # try:
        #     stock_file = parameters["stock_file"]
        # except KeyError:
        #     raise err.ReadConfigError(msg="ERROR! Stock file is not specified in the config file")
        # if stock_file == []:
        #     raise err.ReadConfigError(msg="ERROR! Stock file is not specified in the config file")
        # try:
        #     self.read_stock(stock_file[0])
        # except err.ReadConfigError as e:
        #     raise err.ReadConfigError(msg=e.msg)
        #
        # self.coding_table.recalculate_coding_table()
        #
        # try:
        #     price_file = parameters["prices_file"]
        # except KeyError:
        #     raise err.ReadConfigError(msg="ERROR! Price file is not specified in the config file")
        # if price_file == []:
        #     raise err.ReadConfigError(msg="ERROR! Price file is not specified in the config file")
        # try:
        #     self.read_prices(price_file[0])
        # except err.ReadConfigError as e:
        #     raise err.ReadConfigError(msg=e.msg)

    def read_args(self):
        try:
            job_name = self.args.job_name
            self.parameters_read["job_name"][-2] = job_name
            self.parameters_read["job_name"][2] = True
        except AttributeError:
            pass

        try:
            ncs = self.args.ncs
            self.parameters_read["NCS"][-2] = ncs
            self.parameters_read["NCS"][2] = True
        except AttributeError:
            pass

        try:
            max_block_size = self.args.max_block_size
            self.parameters_read["max_block_size"][-2] = max_block_size
            self.parameters_read["max_block_size"][2] = True
        except AttributeError:
            pass

        try:
            stock = self.args.stock
            self.parameters_read["stock"][-2] = stock
            self.parameters_read["stock"][2] = True
        except AttributeError:
            pass

        try:
            blocks = self.args.blocks
            self.parameters_read["blocks"][-2] = blocks
            self.parameters_read["blocks"][2] = True
        except AttributeError:
            pass
        else:
            self.calculate_blocks = False

        try:
            price_table = self.args.price_table
            self.parameters_read["price_table"][-2] = price_table
            self.parameters_read["price_table"][2] = True
        except AttributeError:
            pass

        if self.args.verbose:
            self.verbose = True

        if self.args.silent:
            self.silent = True

        if self.args.only_blocks:
            self.only_blocks = True

        try:
            self.block_samples = self.args.block_finder_mode[0]
            self.block_min_depth = self.args.block_finder_mode[1]
        except AttributeError:
            pass
        else:
            self.block_finder_mode = True

    def check_parameters(self):
        for key in self.parameters_read:
            if self.parameters_read[key][0]:
                self.parameters_read[key][1] = self.parameters_read[key][-1](self.parameters_read[key][4])
            if self.parameters_read[key][2]:
                self.parameters_read[key][3] = self.parameters_read[key][-1](self.parameters_read[key][5])

    def check_consistency(self):
        if not self.parameters_read['NCS'][1] and not self.parameters_read['NCS'][3]:
            if not self.parameters_read['spectra'][1] or not self.parameters_read['label_types'][1]:
                return False, 'NCS is not specified'

        if self.block_finder_mode:
            if 1 <= self.block_samples <= 9 and 1 <= self.block_min_depth <= 20:
                return True, "OK. Block finder mode."
            else:
                return False, 'Block samples number (1<=s<=9) or minimal depth of search (1<=m<=20) not in specified range'

        if not self.block_finder_mode and not self.only_blocks:
            if not self.parameters_read['price_table'][1] and not self.parameters_read['price_table'][3]:
                return False, 'Price table is not specified or has wrong format'

        if self.calculate_blocks:
            if not self.parameters_read['max_block_size'][1] and not self.parameters_read['max_block_size'][3]:
                return False, 'Max block size is not specified'
        else:
            if not self.parameters_read['blocks'][1] and not self.parameters_read['blocks'][3]:
                return False, 'Elementary blocks file (ELB) is not specified or has wrong format'
        return True, "OK"



        # need ncs if no ncs then read from config
        # mode? if it is find blocks, then we need samples and size
        # if mode is find scheme
        # we need good prices file
        # if stock file specified it should be good
        # then we check if we read or calculate blocks
        # if we read blocks, then we just check if file is good, if it is not specified then we use
        # if amino acids are not specified then we use all 20 as default
        # if calculate blocks then we need max size specified

        pass

    def print_errors(self):
        pass

    @property
    def task_parameters(self):
        task_parameters = {
            "NCS": self.ncs,
            "block_finder_mode": self.block_finder_mode,
            "only_blocks": self.only_blocks,
            "calculate_blocks": self.calculate_blocks,
            "aa_list": self.aa_list,
            "blocks": self.blocks,
            "prices": self.prices,
            "block_samples": self.block_samples,
            "block_min_depth": self.block_min_depth,
            "job_name": self.job_name,
            "max_block_size": self.max_block_size
        }
        return task_parameters

    @property
    def output_parameters(self):
        files = "l"
        if self.calculate_blocks or self.block_finder_mode:
            files += "e"
        if not self.block_finder_mode and not self.only_blocks:
            files += "ps"
        output_parameters = {
            "verbose": self.verbose,
            "silent": self.silent,
            "job_name": self.job_name,
            "files": files
        }
        return output_parameters

    def check_aa_list(self, aa_list):
        for res in aa_list:
            if res not in Constants.RES_TYPES_LIST or res not in Constants.RES_TYPES_THREE:
                return False
        return True

    def read_blocks(self, filename):
        lines = self.read_lines(filename)
        NCS_read = False
        ncs = NC2
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if not NCS_read:
                if line[0] == "NCS":
                    ncs_name = line[1]
                    for each_ncs in LIST_OF_NCS:
                        if each_ncs.name == ncs_name:
                            ncs = each_ncs
                            NCS_read = True
                            break
                i += 1
            else:
                if len(line) == 2:
                    samples = int(line[0])
                    length = int(line[1])
                    patterns = []
                    for j in range(len(length)):
                        patterns.append(lines[i+j+1].strip()[0])
                    block = Scheme("", ncs, samples, patterns)
                    self.add_block(block)
                    i += 1 + length
        return self.blocks

    def add_block(self, block):
        depth_of_scheme = block.size()
        if depth_of_scheme in self.blocks:
            self.blocks[depth_of_scheme].append(block)
        else:
            self.blocks.update({depth_of_scheme: [block]})

    @staticmethod
    def read_lines(filename):
        try:
            with open(filename, 'r', encoding='UTF-8') as f:
                lines = f.readlines()
        except IOError:
            raise err.ReadLinesError(filename)
        new_lines = []
        for line in lines:
            curr_line = line.rstrip()
            if curr_line != "" and curr_line[0] != "#":
                new_lines.append(curr_line.split("#")[0])
        return new_lines