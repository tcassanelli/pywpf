#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
from astropy.time import Time
import argparse

"""
Transform the text file data to the standard numpy binary.
"""

parser = argparse.ArgumentParser()
parser.add_argument(
    "--path-input", "-i",
    help="path to file",
    type=str,
    )

parser.add_argument(
    "--path-output", "-o",
    help="path to file",
    type=str,
    )

args = parser.parse_args()
pi = args.path_input
po = args.path_output

if po is None:
    po = os.path.join(
        os.path.dirname(pi), os.path.basename(pi).split(".txt")[0])

data_mjd = np.genfromtxt(pi, dtype=np.float128)
data_s =  (data_mjd - data_mjd[0]) * 24 * 60 * 60

np.save(po, data_s)
