import logging
from classes.interactive import answer_yes


def create_logger_main(args, default_log_filename):
    logger = logging.getLogger("Logfile_logger")
    logger.setLevel(logging.DEBUG)
    if args.jobname:
        log_filename = "{}.log".format(args.jobname)
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
