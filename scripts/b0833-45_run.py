#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
from astropy import units as u
from astropy.table import Table
import matplotlib.pyplot as plt
import argparse
import pywpf

"""
Simple script to run PyWPF on Crab pulsar data.
The script will run and generate the output data in the working directory.
"""

parser = argparse.ArgumentParser()
parser.add_argument(
    "--path-input", "-i",
    help="path to file",
    type=str,
    )
args = parser.parse_args()
times_path = args.path_input

# B0833-45
dt = 0.002793
T_init = 0.089367
iteration = 100
delta = 10e-9

work_dir = os.path.dirname(times_path)
basename = os.path.basename(times_path)
name = os.path.splitext(basename)[0]

M = 20
pywpf.pca_folding(
    times_path=times_path,
    dt=dt,
    T_init=T_init,
    iteration=iteration,
    delta=delta,
    num_div=[M],
    merit_func=pywpf.merit1,
    region_order=3,
    work_dir=work_dir
    )

# beyond this point we take data and plot it
n_counts = np.load(times_path).size
pywpf_path = os.path.join(work_dir, "pywpf_out")
data = np.load(os.path.join(pywpf_path, f"{name}-000", f"M{M}.npz"))
pywpf_info = os.path.join(pywpf_path, f"{name}-000", 'info.dat')
info = Table.read(pywpf_info, format="ascii")

# Change in the future, all of the merits have same time sample
factor = 89.3 * u.ms
time_x = np.linspace(
    info['T_init'][0] - info['delta'][0] * info['iter'][0] / 2,
    info['T_init'][0] + info['delta'][0] * info['iter'][0] / 2,
    info['iter'][0],
    endpoint=False
    ) * u.s - factor

ev = data["EVALW"][:, 0]
s = data["SW"][:, 0]
merit = data["MERIT"]

fig, ax = plt.subplots()

ax.plot(
    time_x.to(u.us), merit,
    label="".join(("$\\xi_\\ell(M=", f"{M}", ")$")), c="k", alpha=0.5, lw=1,
    )
ax.plot(
    time_x.to(u.us), s,
    label="$\\left|s_{\\ell, 1}\\right|$", c="k", alpha=1, ls='--')
ax.plot(
    time_x.to(u.us), ev,
    label="$\\lambda_{\\ell, 1}$", c="k", alpha=1, ls=":")

ax.set_xlim(time_x.min().to_value(u.us), time_x.max().to_value(u.us))
ax.grid(True)
ax.legend(loc="upper right")
ax.set_xlabel(
    "".join((
        f"Time ms + ${factor.to_value(u.ms)}$ ms",
        f" (Total counts: {n_counts}; $M={M}$)"
        ))
    )
ax.set_ylabel("Eigenvalue, scalar and merit function")
fig.tight_layout()
fig.savefig(os.path.join(pywpf_path, f"{name}_run.pdf"))
