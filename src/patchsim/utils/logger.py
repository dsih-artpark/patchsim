import logging
import os

def setup_logger(log_file="output/logs/io.log"):
    """
    Set up a logger to log messages to a file and the console.

    Args:
        log_file (str): Path to the log file.
    """
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    # Configure the logger
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger("PatchSimLogger")