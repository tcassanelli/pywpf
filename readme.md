# PCA Folding Analysis
### Author: Tomas Cassanelli
Staring date: 16-Sep-2016

## Summary
The purpose of this program is to see the effectivity and efficiency that a Principal Component Analysis (PCA) folding analysis has compared to the **Epoch Folding**, which is the standard used for pulsar period computations, with observations in the optical regime.

## Program
The program is subdivided in two scripts, the first one, `pca_analysis.py`, it has all the functions and calculations and the second one, `pca_run.py`, runs the program and stores the generated matrices or plots in the working directory. Another extention is meant to be build which will be a GUI (soon).

### Inputs
1. CSV file which contains the time array that comes from the instrument. It is the time array that comes straight from the astronomical instrument when a pulsar is observed. Basically one column vector.
2. The bintime selected, in the program is also called `dt`.
3. The `period_start`, or the given period to start the iteration. Depending on how noisy or clean is the time array this can be computed or approximated since the beginning using a FFT. If the signal is not clear it is better to start with a close value to the real period. Notice that since these periods are well known, starting close to the real value is not a problem.
4. Number of iterations, `iter1`, `iter2`, there are two big iterations loops, the second one should be the one with more iteration steps.
5. `delta1`, `delta2`, are the amount of time added or substracted depending on the iterations. They are of the order of the ns.
6. The `num_div`, which is in fact the number of how the time array is going to be divided. Also corresponds to the number of rows in the waterfall diagram and to the eigenvalues from the PCA analysis. Please refer to [A Tutorial on Principal Component Analysis](https://arxiv.org/abs/1404.1100) by Jonathon Shlens.

### Method
To obtain the best period selection, the PCA method is applied. For this, a folding is done to the waterfall diagram, then the covariance matrix is computed, diagonalized and normalized with mean zero and standard deviation 1. The output gives the values for the Principal Component (PC) vectors, **P**, and the eigenvalues, **V**. Depending on the number of rows in the waterfall will be obtained different number of dimensions. The number of dimension it is also number of eigenvalues. For one iteration the result would be,

<p align="center">
<img src="https://github.com/tcassanelli/PCA-Folding/blob/master/images/eq1.png" alt="Eigenvector and eigenvalues" width="300">
</p>

Where M is the number of dimentions or `num_div`. Another important quantity to find the best period corresponds to the scalar, S. The scalar corresponds to the dot product between the eigenvectors and the hyperdiagonal unitary vector. This is

![Hyperdiagonal unitary vector](https://github.com/tcassanelli/PCA-Folding/blob/master/images/eq2.png)

Adding for N iterations (`iter1` or `iter2`) to a single matrix for the case of the scalar [**S**] and the eigenvalue [**V**],

![Scalar and eigenvector](https://github.com/tcassanelli/PCA-Folding/blob/master/images/eq3.png)








### Outputs 

The principal outputs of the program will be the second and first periods from the iterations. It is possible to obtain other values from the residual data.

1. `period_final1`, `period_final2`, correspond to the best period found using the method. In an ascii file.
2. Eigenvalues, in a `.npy` output format.
3. Scalars (hyperdiagonal unitary vector times the eigenvector), in a `.npy` output format.
4. MSTEV (Maximum Scalar Times EigenValue), which is the selected merit function for the analysis. Same output as before, `.npy`.
