# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os
import ruamel.yaml as yaml

"""
This __init__.py file is mainly used to set up the package wide config.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "0.2"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "Apr 2018"


SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".quotas_config.yaml")


def _load_pmg_settings():
    try:
        with open(SETTINGS_FILE, "rt") as f:
            d = yaml.safe_load(f)
    except FileNotFoundError:
        d = {}

    return d

SETTINGS = _load_pmg_settings()