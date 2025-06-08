import logging


def lr_multiply_ratio(initial_lr, lr_ratio):
    def lr_multiplier(idx):
        return initial_lr * lr_ratio**idx

    return lr_multiplier


def create_logger(name: str) -> logging.Logger:
    """
    Creates a logger with a stream handler and two file handlers.
    The stream handler prints to the screen depending on the value of `quiet`.
    One file handler (verbose.log) saves all logs, the other (quiet.log) only saves important info.
    :param save_dir: The directory in which to save the logs.
    :return: The logger.
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.propagate = False

    # create file handler which logs even debug messages
    fh = logging.FileHandler(f"{name}/output.log")
    fh.setLevel(logging.DEBUG)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    # create formatter and add it to the handlers
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger
