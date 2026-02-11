import logging
import cmdstanpy

def setup_loggers():
    birdman_logger = logging.getLogger("birdman")
    birdman_logger.setLevel(logging.INFO)
    sh = logging.StreamHandler()
    formatter = logging.Formatter(
        "[%(asctime)s - %(name)s - %(levelname)s] ::  %(message)s"
    )
    sh.setFormatter(formatter)
    sh.setLevel(logging.INFO)
    birdman_logger.addHandler(sh)

    cmdstanpy_logger = cmdstanpy.utils.get_logger()
    cmdstanpy_logger.setLevel(logging.INFO)
    cmdstanpy_logger.addHandler(sh)
    for h in cmdstanpy_logger.handlers:
        h.setFormatter(formatter)

    return birdman_logger
