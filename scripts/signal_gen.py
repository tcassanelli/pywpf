import os
import numpy as np

# generating plain and simple signal, but how? things you need to consider:
# number of total data points
# the max and min delta between those data points, pico seconds?
# Initial and final time
# how many pulses shoud I introduce in this signal, same as the vela pulse?


def noise_to_signal(original_signal, level):
    """
    level is how many data points of noise will be added to the signal
    original_signal is basically the time array.
    """

    x_size = int(original_signal.size * level)

    X = np.random.random_sample(x_size)

    ti = original_signal[0]
    tf = original_signal[-1]

    time_noise = X * (tf - ti) + ti

    new_signal = np.hstack((original_signal, time_noise))

    new_signal_sorted = np.sort(new_signal)

    return new_signal_sorted


if __name__ == '__main__':

    # Now generating the data files
    # 5550238 original size from Enrico
    data = np.load('../../data_B0833-45/data_B0833-45.npy')[:5550238]

    level_list = np.array([0, 25, 50, 75, 100, 200, 250, 275, 300, 325]) / 100
    for k, level in enumerate(level_list):
        noise_data = noise_to_signal(original_signal=data, level=level)
        np.save(
            file=os.path.join("../../data_B0833-45", f"n{int(level * 100)}"),
            arr=noise_data
            )
        del(noise_data)
        print(f'{k} data set completed')
