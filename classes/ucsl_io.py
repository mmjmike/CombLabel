from cl_errors import errors as err
from classes.constants import Constants, NCS, ELB
# from classes.search_objects import ELB
import re
import os
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
        self.products_filename = self.job_name + "_products.txt"
        self.opened_files = {}
        self.open_files()
        self.first_products = True
        self.first_bl_stats = True

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

    def write_data(self, output, files="l", timer=1):
        types = list(files)
        if "l" in types and "l" in self.opened_files:
            if timer :
              timestamp = "\n" + datetime.datetime.now().strftime('[Date,Time = %Y-%m-%d_%H:%M:%S]')
            else :
              timestamp = ""
            self.opened_files["l"].write(timestamp + "\n" + output)
            self.opened_files["l"].flush()
        if "e" in types and "e" in self.opened_files:
            self.opened_files["e"].write(output)
            self.opened_files["e"].flush()
        if not self.silent:
            if "c" in types or self.verbose:
                print(output)

    def close_files(self):
        for key in self.opened_files:
            self.opened_files[key].close()

    def write_best_scheme(self, best_scheme):
        output = '________________________\nBest scheme\n[NCS = {}]\n'.format(best_scheme.scheme.ncs.name)
        output += "[solution]\nRes"
        for i in range(best_scheme.samples):
            output += ",S{}".format(i + 1)
        output += "\n"
        for res in best_scheme.residues:
            output += res + ", " + ", ".join(list(best_scheme.label_dict[res])) + "\n"
        output += "[price]\n{}\n\n".format(best_scheme.price)
        output += "Blocks used for this scheme:\n"
        for block in best_scheme.blocks:
            output += "[ELB samples = {} patterns = {}]\n".format(block.samples, len(block.patterns)) \
                      + str(block)
        with open(self.job_name + "_best_scheme.txt", mode="w") as f:
            f.write(output)
        self.write_data(output, files="cl")

    def write_product_stats(self, stats):
        output = ""
        for samples in sorted(stats):
            output += "-------------\n"
            output += "Products in {} samples:\n".format(samples)
            for product_type in sorted(stats[samples]):
                output += product_type + "\n"
                for status in sorted(stats[samples][product_type]):
                    output += "{:4} scheme(s): {}\n".format(stats[samples][product_type][status], status)
        with open(self.job_name + "_product_stats.txt", mode="w") as f:
            f.write(output)
        self.write_data(output, files="cl")

    def write_products(self, products, samples):
        if not products:
            return
        product_schemes = 0
        output = ""
        if self.first_products:
            output += "[NCS = {}]\n".format(products[0].ncs.name)
        output += self.make_block_stats(products[0].blocks)
        output += "\n-----------------------\nProducts calculated for {} samples:\n".format(samples)
        for product in products:
            output += str(product) + ": {} scheme(s)\n".format(len(product))
            product_schemes += len(product)
        output += "{} product types calculated\n".format(len(products))
        output += "{} total labeling schemes to check\n".format(product_schemes)

        mode = "a"
        if self.first_products:
            mode = "w"
            self.first_products = False
        with open(self.job_name + "_products.txt", mode=mode) as f:
            f.write(output)
        self.write_data(output, files="cl")

    def write_block_stats(self, blocks):
        output = make_block_stats(blocks)
        mode = "a"
        if self.first_bl_stats:
            mode = "w"
            self.first_bl_stats = False
        with open(self.job_name + "_elb_stats.txt", mode=mode) as f:
            f.write(output)
        self.write_data(output, files="cl")


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
        except TypeError:
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
        try:
            return len(job_name) > 0
        except TypeError:
            return False

    def check_ncs(self, ncs):
        try:
            return ncs.upper() in Constants.NCS_NAMES
        except AttributeError:
            return False

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
        except (FileNotFoundError, TypeError):
            return False
        return self.read_prices(prices_file)[0]

    def read_prices(self, prices_file):
        print("reading prices from file {}".format(prices_file))
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
            if residue_type in Constants.RES_TYPES_THREE:
                residue_type = Constants.TO_ONE_LETTER_CODE[residue_type]
            if residue_type not in Constants.RES_TYPES_LIST:
                msg = "Wrong residue type '{}' in prices file '{}'".format(residue_type, prices_file)
                return [False, msg, prices]
            prices.update({residue_type: curr_dict})
        print(prices)
        return [True, msg, prices]

    def check_stock_file(self, stock_file):
        try:
            f = open(stock_file, "r")
            f.close()
        except (FileNotFoundError, TypeError):
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
            print(residue_type)
            if residue_type in Constants.RES_TYPES_THREE:
                residue_type1 = Constants.TO_ONE_LETTER_CODE[residue_type]
                print("{} -> {}".format(residue_type,residue_type1))
                residue_type = residue_type1
            if residue_type not in Constants.RES_TYPES_LIST:
                msg = "Wrong residue type '{}' in stock file '{}'".format(residue_type, stock_file)
                return [False, msg, label_dict]
            label_dict.update({residue_type: curr_list})
        return [True, msg, label_dict]

    def check_block_size(self, block_size):
        try:
            size = int(block_size)
        except (ValueError, TypeError):
            return False
        return 0 < size < 10

    def check_blocks_file(self, block_file):
        try:
            f = open(block_file, "r")
            f.close()
        except (FileNotFoundError, TypeError):
            return False
        result = self.read_block_file(block_file)
        return result[0]

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
                patterns = lines[i+1:i+1+patterns_num]
                good_block = True
                for pattern in patterns:
                    if len(pattern) != samples_num:
                        msg = "Warning! The number of samples in block {} doesn't match " \
                              "\nwith specified value in blocks file '{}'.".format(blocks_num, block_file)
                        good_block = False
                    for label_type in pattern:
                        if label_type not in Constants.TYPES:
                            msg = "Warning! Unknown labeling type '{}' used in block {}" \
                                  "\n in blocks file '{}'.".format(label_type, blocks_num, block_file)
                            good_block = False
                            break
                if good_block:
                    block = Scheme("", ncs, samples_num, patterns)
                    if samples_num not in blocks:
                        blocks.update({samples_num: {patterns_num: [block]}})
                    else:
                        if patterns_num not in blocks[samples_num]:
                            blocks[samples_num].update({patterns_num: [block]})
                        else:
                            blocks[samples_num][patterns_num].append(block)
                i += 1 + patterns_num
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
            print(msg)

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
                self.blocks = self.read_block_file(blocks_file)[2]
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
                self.parameters_read[key][4] = parameters[key][0]

        if self.parameters_read["blocks"][0]:
            self.calculate_blocks = False

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
            if blocks:
                self.parameters_read["blocks"][-2] = blocks
                self.parameters_read["blocks"][2] = True
                self.calculate_blocks = False
        except AttributeError:
            pass

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
        except TypeError:
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
            print(self.parameters_read['blocks'])
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
            files += "p"
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
        ncs = Constants.NC2
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if not NCS_read:
                if line[0] == "NCS":
                    ncs_name = line[1]
                    for each_ncs in Constants.LIST_OF_NCS:
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


def make_block_stats(blocks):
    output = "Blocks used:\n"
    total_blocks = 0
    samples = list(blocks.keys())
    samples.sort()
    blocks_by_samples = {}
    for samples_n in samples:
        pattern_numbers = list(blocks[samples_n].keys())
        pattern_numbers.sort()
        blocks_in_samples = 0
        for patterns_n in pattern_numbers:
            blocks_n = len(blocks[samples_n][patterns_n])
            total_blocks += blocks_n
            blocks_in_samples += blocks_n
            output += "{:2} X  {:2} - {} block(s)\n".format(samples_n, patterns_n, blocks_n)
        blocks_by_samples.update({samples_n: blocks_in_samples})
    output += "-----------\n"
    for samples_n in samples:
        output += "{:2} X  . - {} block(s)\n".format(samples_n, blocks_by_samples[samples_n])
    output += " TOTAL - {} blocks\n".format(total_blocks)
    output += "-----------\n"
    return output


def read_blocks(block_file, logger=None):
    lines = read_lines(block_file)
    i = 0
    deuterated = False
    deuterated_found = False
    result = ''
    blocks = {}
    blocks_num = 0
    NCS_found = False
    ncs_name = ''
    ncs_regular = Constants.ncs_re
    while i < len(lines):
        if not NCS_found:
            ncs_match = ncs_regular.match(lines[i])
            if ncs_match:
                ncs_name = ncs_match.group(1).upper()
                NCS_found = True
                i += 1
                continue
        if not deuterated_found:
            deuterated_match = Constants.deuterated_re.match(lines[i])
            if deuterated_match:
                deuteration = deuterated_match.group(1).upper()
                deuterated_found, deuterated = extract_deuterated(deuteration)
                i += 1
                continue
        if NCS_found:
            elb_match = Constants.elb_re.match(lines[i])
            if elb_match:
                samples_num = int(elb_match.group(1))
                patterns_num = int(elb_match.group(2))
                blocks_num += 1
                patterns = lines[i + 1:i + 1 + patterns_num]
                good_block = True
                for pattern in patterns:
                    if not good_block:
                        break
                    if len(pattern) != samples_num:
                        result = "Warning! The number of samples in block {} doesn't match " \
                              "\nwith specified value in blocks file '{}'.".format(blocks_num, block_file)
                        print(result)
                        if logger:
                            logger.warning(result)
                        good_block = False
                        break
                    for label_type in pattern:
                        if label_type not in Constants.TYPES:
                            result = "Warning! Unknown labeling type '{}' used in block {}" \
                                  "\n in blocks file '{}'.".format(label_type, blocks_num, block_file)
                            if logger:
                                logger.warning(result)
                            else:
                                print(result)
                            good_block = False
                            break
                if good_block:
                    block = ELB(patterns, ncs_name, deuterated)
                    if samples_num not in blocks:
                        blocks.update({samples_num: {patterns_num: [block]}})
                    else:
                        if patterns_num not in blocks[samples_num]:
                            blocks[samples_num].update({patterns_num: [block]})
                        else:
                            blocks[samples_num][patterns_num].append(block)
                i += 1 + patterns_num
        else:
            i += 1
    return result, blocks, ncs_name, deuterated


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


def extract_labels(line):
    extracted_labels = [label.rstrip().upper() for label in line.split(",")]
    for label in extracted_labels:
        if label not in Constants.TYPES_NAMES:
            labels_names = ", ".join(Constants.TYPES_NAMES)
            print("Error! Label type '{}' is not in the following list:\n{}".format(label, labels_names))
            return []
    labels = []
    for label in Constants.BASIC_TYPES:
        if label.name in extracted_labels:
            labels.append(label)
    return labels


def extract_spectra(line):
    extracted_spectra = [spectrum.rstrip() for spectrum in line.split(",")]
    for spectrum in extracted_spectra:
        if spectrum not in Constants.SPECTRA_NAMES:
            spectra_names = ", ".join(Constants.SPECTRA_NAMES)
            print("Error! Spectrum '{}' is not in the following list:\n{}".format(spectrum, spectra_names))
            return []
    spectra = []
    for spectrum in Constants.basic_spectra:
        if spectrum.name in extracted_spectra:
            spectra.append(spectrum)
    return spectra


def extract_deuterated(deuteration):
    result = False
    deuterated = False
    d = deuteration.upper()
    if d == "TRUE" or d == "YES" or d == "Y" or d == "1":
        result = True
        deuterated = True
    elif d == "FALSE" or d == "NO" or d == "N" or d == "0":
        result = True
    return result, deuterated


def read_ncs(filename):
    lines = read_lines(filename)
    ncs_name_read = False
    labels_read = False
    spectra_read = False
    deuterated_read = False
    ncs_name = ''
    spectra = []
    labels = []
    deuterated = False


    for line in lines:
        if not ncs_name_read:
            ncs_name_search = ncs_re.search(line)
            if ncs_name_search is not None:
                ncs_name = ncs_name_search.group()[0]
                if ncs_name:
                    ncs_name_read = True
                    continue
        if not labels_read:
            labels_search = labels_re.search(line)
            if labels_search is not None:
                labels_line = labels_search.group()[0]
                labels = extract_labels(labels_line)
                if labels:
                    labels_read = True
                    continue
        if not spectra_read:
            spectra_search = spectra_re.search(line)
            if spectra_search is not None:
                spectra_line = spectra_search.group()[0]
                spectra = extract_spectra(spectra_line)
                if spectra:
                    spectra_read = True
                    continue
        if not deuterated_read:
            deuterated_search = deuterated_re.search(line)
            if deuterated_search is not None:
                deuterated_line = deuterated_search.group()[0]
                deuterated = extract_deuterated(deuterated_line)
                deuterated_read = True
                continue

    if labels_read and ncs_name_read and spectra_read:
        return NCS(ncs_name, labels, spectra, deuterated)
    else:
        return None


def add_dir(path_list, new_dir):
    new_path_list = []
    for path in path_list:
        new_path_list.append(path)
        new_path_list.append(os.path.join(path, new_dir))
    return new_path_list


def add_file(path_list, filename):
    ncs_ext = ".ncs"
    new_path_list = []
    for path in path_list:
        new_path_list.append(os.path.join(path, filename))
        new_path_list.append(os.path.join(path, filename + ncs_ext))
    return new_path_list


def check_ncs_path(full_path):
    # print("Checking path: '{}'".format(full_path))
    if os.path.isfile(full_path):
        # print("File exists!")
        return True
    else:
        # print("File not found :(")
        return False


def add_to_file_name(filename, addition):
    path, only_name = os.path.split(filename)
    name, ext = os.path.splitext(only_name)
    name += addition
    only_name = name + ext
    return os.path.join(path, only_name)


def find_good_ncs_paths(ncs_name, script_path=os.path.split(os.path.realpath(__file__))[0]):
    ncs_dir = "NCS"
    paths = [os.getcwd(), script_path]
    paths = add_file(add_dir(paths, ncs_dir), ncs_name)
    good_paths = []
    for path in paths:
        if check_ncs_path(path):
            good_paths.append(path)
    return good_paths


def write_blocks(blocks, ncs_name, filename, deuterated):
    output = "[NCS = {}]\n".format(ncs_name)
    output+= "[Deuterated = {}]\n".format(deuterated)
    blocks_samples = list(blocks.keys())
    blocks_samples.sort()
    for samples_num in blocks_samples:
        patterns_numbers = list(blocks[samples_num].keys())
        patterns_numbers.sort()
        for patterns_num in patterns_numbers:
            for block in blocks[samples_num][patterns_num]:
                output += "[ELB samples = {} patterns = {}]\n".format(block.samples, len(block.patterns)) \
                          + str(block)
    output+="\n"
    with open(filename, mode="w") as f:
        f.write(output)


def find_ncs(ncs_name, script_path):
    paths = find_good_ncs_paths(ncs_name, script_path)
    for path in paths:
        ncs = read_ncs_file(path)
        if ncs:
            msg = "NCS '{}' read from '{}'".format(ncs_name, path)
            return ncs, msg
    msg = "Error! NCS file not found"
    return None, msg


def read_ncs_file(filename):
    name = ''
    deuterated = False
    deuterated_read = False
    labels = []
    spectra = []
    with open(filename, "r") as f:
        for line in f:
            if not name:
                name_match = Constants.ncs_re.match(line)
                if name_match:
                    name = name_match.group(1)
                    continue
            if not labels:
                labels_match = Constants.labels_re.match(line)
                if labels_match:
                    labels = extract_labels(labels_match.group(1))
                    continue
            if not spectra:
                spectra_match = Constants.spectra_re.match(line)
                if spectra_match:
                    spectra = extract_spectra(spectra_match.group(1))
                    continue
            if not deuterated_read:
                deuterated_match = Constants.ncs_re.match(line)
                if deuterated_match:
                    deuteration = deuterated_match.group(1).upper()
                    deuterated_found, deuterated = extract_deuterated(deuteration)
                    continue
    if name and labels and spectra:
        return NCS(name, spectra, labels, deuterated)
    return None


def write_ncs_stamp(ncs):
    return "[NCS = {}]\n[Deuterated = {}]".format(ncs.name, ncs.deuterated)


def write_best_scheme(best_scheme, filename):
    output = '________________________\nBest scheme\n[NCS = {}]\n'.format(best_scheme.scheme.ncs_name)
    output += "[solution]\nRes"
    for i in range(best_scheme.samples):
        output += ",S{}".format(i + 1)
    output += "\n"
    for res in best_scheme.residues:
        output += res + ", " + ", ".join(list(best_scheme.label_dict[res])) + "\n"
    output += "[price]\n{}\n\n".format(best_scheme.price)
    output += "Blocks used for this scheme:\n"
    for block in best_scheme.blocks:
        output += block.full_str() + "\n"
    with open(filename, mode="w") as f:
        f.write(output)
    return output


def write_product_stats(stats, filename):
    output = ""
    for samples in sorted(stats):
        output += "-------------\n"
        output += "Products in {} samples:\n".format(samples)
        for product_type in sorted(stats[samples]):
            output += product_type + "\n"
            for status in sorted(stats[samples][product_type]):
                output += "{:4} scheme(s): {}\n".format(stats[samples][product_type][status], status)
    with open(filename, mode="w") as f:
        f.write(output)
    return output


def write_products(products, samples, filename, mode='w'):
    if not products:
        return
    product_schemes = 0
    output = ""
    # if mode == 'w':
    #     output += write_ncs_stamp(products[0].ncs)
    output += make_block_stats(products[0].blocks)
    output += "\n-----------------------\nProducts calculated for {} samples:\n".format(samples)
    for product in products:
        output += str(product) + ": {} scheme(s)\n".format(len(product))
        product_schemes += len(product)
    output += "{} product types calculated\n".format(len(products))
    output += "{} total labeling schemes to check\n".format(product_schemes)
    with open(filename, mode=mode) as f:
        f.write(output)
    return output


def read_prices(prices_file):
    msg = ""
    deuterated_found = False
    prices = {}
    lines = read_lines(prices_file)
    if len(lines) < 3:
        msg = "Prices file '{}' is too short".format(prices_file)
        return prices, msg
    d = []
    line_length = 0
    first_line = True
    for line in lines:
        deuterated_match = Constants.deuterated_re.match(line)
        if deuterated_match:
            deuteration = deuterated_match.group(1).upper()
            deuterated_found, deuterated = extract_deuterated(deuteration)
            prices.update({"Deuterated" : deuterated})
            continue
        else:
          s = [x.strip() for x in line.split(",")]
          if first_line:
              line_length = len(s)
              first_line = False
          else:
              if len(s) != line_length:
                  msg = "Not equal length of lines in prices file '{}'".format(prices_file)
                  return prices, msg
          d.append(s)
    price_label_types = d[0][1:]
    for l_type in price_label_types:
        if l_type.upper() not in Constants.TYPES_NAMES:
            msg = "Incorrect label type {} in prices file '{}'".format(l_type, prices_file)
            return prices, msg
    if len(price_label_types) < 2:
        msg = "Too few labeling types specified in prices file '{}'".format(prices_file)
        return prices, msg

    for i in range(len(d) - 1):
        curr_dict = {}
        for j in range(len(price_label_types)):
            try:
                price = float(d[i + 1][j + 1])
            except ValueError:
                msg = "ERROR! Price must be set in digits"
                msg += "\nPlease check price file '{}' (row {}; col {})".format(prices_file,
                                                                                i + 2, j + 2)
                return {}, msg
            curr_dict.update({price_label_types[j]: price})
        residue_type = d[i + 1][0]
        if residue_type not in Constants.RES_TYPES_LIST and residue_type not in Constants.RES_TYPES_THREE:
            msg = "Wrong residue type '{}' in prices file '{}'".format(residue_type, prices_file)
            return {}, msg
        if residue_type in Constants.RES_TYPES_THREE:
            residue_type1 = Constants.TO_ONE_LETTER_CODE[residue_type]
            residue_type = residue_type1
        prices.update({residue_type: curr_dict})
    return prices, msg
