
import numpy as np

# generating plain and simple signal, but how? things you need to consider:
# number of total data points
# the max and min delta between those data points, pico seconds?
# Initial and final time
# how many pulses shoud I introduce in this signal, same as the vela pulse?


def noise_to_signal(old_signal, level):
    """
    level is how many data points of noise will be added to the signal
    old_signal is basically the time array.
    """

    x_size = int(old_signal.size * level)

    X = np.random.random_sample(x_size)

    ti = old_signal[0]
    tf = old_signal[-1]

    time_noise = X * (tf - ti) + ti

    new_signal = np.hstack((old_signal, time_noise))

    new_signal_sorted = np.sort(new_signal)

    return new_signal_sorted


# Now generating the data files
data = np.load('../../data_pulsar/n0.npy')

for l in np.arange(0, 1.05, 0.05):

    noise_data = noise_to_signal(old_signal=data, level=l)

    np.save(
        file='n' + str(int(l * 100)),
        arr=noise_data
        )
