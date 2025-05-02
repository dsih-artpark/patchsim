import logging
import os
import platform
import sys
from datetime import datetime

def setup_logger(model_name, config, num_patches, patches, base_model):
    """
    Set up a logger to log messages to a file and the console, and log system/run details.
    Args:
        model_name (str): Name of the model.
        config (dict): Configuration dictionary.
        num_patches (int): Number of patches.
        patches (list): List of patch names.
        base_model (CompartmentalModel): Base model object.
    Returns:
        logging.Logger: Configured logger.
    """
    log_dir = os.path.join(config['OutputDir'], "logs")
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f"{model_name}_run_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
    logger = logging.getLogger("PatchSimLogger")
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    fh = logging.FileHandler(log_file)
    fh.setFormatter(formatter)
    logger.handlers = []
    logger.addHandler(fh)
    # Log system and run details
    logger.info(f"Model: {model_name}")
    logger.info(f"Python version: {sys.version}")
    logger.info(f"Platform: {platform.platform()}")
    logger.info(f"Parameters: {base_model.parameters}")
    # Parameter agnostic positivity check
    for param, value in base_model.parameters.items():
        try:
            if float(value) <= 0:
                logger.warning(f"Parameter '{param}' has non-positive value: {value}")
        except Exception:
            logger.warning(f"Parameter '{param}' could not be checked for positivity (value: {value})")
    logger.info(f"PatchFile: {config['PatchFile']}")
    logger.info(f"SeedFile: {config['SeedFile']}")
    logger.info(f"NetworkFile: {config['NetworkFile']}")
    logger.info(f"OutputDir: {config['OutputDir']}")
    logger.info(f"TMax: {config['TMax']}")
    logger.info(f"Num patches: {num_patches}")
    logger.info(f"Patch list: {patches}")
    logger.info(f"Base model: compartments={base_model.compartments}, transitions={base_model.transitions}, parameters={base_model.parameters}")
    return logger