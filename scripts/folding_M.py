
import os
import numpy as np
import matplotlib.pyplot as plt
from stingray.pulse.pulsar import fold_events
from stingray.pulse import epoch_folding_search
from astropy import units as u

pypcaf_path_crab = "/scratch/v/vanderli/cassane/data_B0833-45"
times = np.load(os.path.join(pypcaf_path_crab, 'data_B0833-45.npy'))

delta = 10e-9
T_init = 0.08936709
dt = 0.0002793
iteration = 1000

number_divisons = 20
N = round(T_init / dt)
# bins = np.linspace(0, 1, N + 1)

T_iteration = np.linspace(
    start=T_init - delta * iteration / 2,
    stop=T_init + delta * iteration / 2,
    num=iteration,
    endpoint=False,
    dtype=times.dtype
    )

efstat = epoch_folding_search(times, 1 / T_iteration, nbin=N)[1]

times_M = np.array_split(times, number_divisons)
del(times)
efstat_M = np.zeros(
    shape=(T_iteration.size, number_divisons),
    dtype=times_M[0].dtype
    )
for M in range(number_divisons):
    efstat_M[:, M] = epoch_folding_search(
        times_M[M], 1 / T_iteration, nbin=N
        )[1]

t_factor = 89.3 * u.ms
T_iteration = T_iteration * u.s

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot((T_iteration - t_factor).to_value(u.us), efstat / efstat.max(), lw=.5)
ax.plot(
    (T_iteration - t_factor).to_value(u.us),
    efstat_M.mean(axis=1) / efstat.max(),
    lw=.5
    )
ax.set_xlabel(f"Time [us] + {t_factor}")
fig.savefig(os.path.join(pypcaf_path_crab, "comparison_B0833.pdf"))
plt.close(fig)


# running the crab

pypcaf_path_crab = "/scratch/v/vanderli/cassane/data_B0531+21"
times = np.load(os.path.join(pypcaf_path_crab, 'data_B0531+21.npy'))

dt = 0.0001
T_init = 0.0336372543236884
iteration = 10000
delta = 1e-9

number_divisons = 20
N = round(T_init / dt)

T_iteration = np.linspace(
    start=T_init - delta * iteration / 2,
    stop=T_init + delta * iteration / 2,
    num=iteration,
    endpoint=False,
    dtype=times.dtype
    )

efstat = epoch_folding_search(times, 1 / T_iteration, nbin=N)[1]

times_M = np.array_split(times, number_divisons)
del(times)
efstat_M = np.zeros(
    shape=(T_iteration.size, number_divisons), 
    dtype=times_M[0].dtype
    )
for M in range(number_divisons):
    efstat_M[:, M] = epoch_folding_search(
        times_M[M], 1 / T_iteration, nbin=N
        )[1]

t_factor = 89.3 * u.ms
T_iteration = T_iteration * u.s

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot((T_iteration - t_factor).to_value(u.us), efstat / efstat.max(), lw=.5)
ax.plot(
    (T_iteration - t_factor).to_value(u.us),
    efstat_M.mean(axis=1) / efstat.max(),
    lw=.5
    )
ax.set_xlabel(f"Time [us] + {t_factor}")
fig.savefig(os.path.join(pypcaf_path_crab, "comparison_B0531.pdf"))
plt.close(fig)
