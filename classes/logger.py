import logging
import os
from classes.interactive import answer_yes


def create_logger_main(args, default_log_filename):
    logger = logging.getLogger("Logfile_logger")
    logger.setLevel(logging.DEBUG)
    if args.jobname:
        log_filename = "{}.log".format(args.jobname)
    else:
        try:
            ncs_name = args.ncs.split(".")[0]
        except AttributeError:
            ncs_name = ""
        try:
            seq_name = args.ncs.split(".")[0]
        except AttributeError:
            seq_name = ""
        if seq_name and ncs_name:
            log_filename = "{}_{}.log".format(seq_name, ncs_name)
        else:
            log_filename = default_log_filename
        if os.path.isfile(log_filename):
            input_key = input("Logfile '{}' (default name) already exists."
                              "You can either overwrite this file, or specify the jobname by using "
                              "--jobname (-j) argument. Do you want to overwrite? (Y[es]/N[o]) > ".format(log_filename))
            if not answer_yes(input_key):
                exit()
    logfile_handler = logging.FileHandler(log_filename, mode='w')
    logfile_format = logging.Formatter('%(asctime)s:%(message)s')
    logfile_handler.setFormatter(logfile_format)
    logfile_handler.setLevel(logging.DEBUG)
    stream_handler = logging.StreamHandler()
    stream_level = logging.INFO
    if args.verbose:
        stream_level = logging.DEBUG
    if args.silent:
        stream_level = logging.CRITICAL
    stream_handler.setLevel(stream_level)
    logger.addHandler(logfile_handler)
    logger.addHandler(stream_handler)
    return logger
