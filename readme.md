# PCA Folding Analysis
### Author: Tomas Cassanelli
Staring date: 16-Sep-2016

## Summary
The purpose of this program is to see the effectivity and efficiency that a PCA folding analysis has compared to the **Epoch Folding**, which is the standard used for pulsar period computations, with observations in the optical regime.

## Program
The program is subdivided in two scripts, the first one, `pca_analysis.py`, it has all the functions and calculations and the second one, `pca_run.py`, runs the program and stores the generated matrices or plots in the working directory. Another extention is meant to be build which will be a GUI (soon).

### Inputs
1. CSV file which contains the time array that comes from the instrument. It is the time array that comes straight from the astronomical instrument when a pulsar is observed.
2. The bintime selected, in the program is also called `dt`.
3. The `period_start`, or the given period to start the iteration. Depending on how noisy or clean is the time array this can be computed or approximated since the beginning using a FFT. If the signal is not clear it is better to start with a close value to the real period. Notice that since these periods are well known, starting close to the real value is not a problem.
4. Number of iterations, `iter1`, `iter2`, there are two big iterations loops, the second one should be the one with more iteration steps.
5. `delta1`, `delta2`, are the amount of time added or substracted depending on the iterations. They are of the order of the ns.
6. The `num_div`, which is in fact the number of how the time array is going to be divided. Also corresponds to the number of rows in the waterfall diagram and to the eigenvalues from the PCA analysis. Please refer to [A Tutorial on Principal Component Analysis](https://arxiv.org/abs/1404.1100) by Jonathon Shlens.

### Method

To obtain the best period selection, the PCA method is applied. For this, a folding is applied to the waterfall diagram, then the covariance matrix is computated and diagonalized. The output gives the values for the PC vectors and V, the eigenvalues. Depending on the number of rows in the waterfall will be obtained different number of dimensions. This will be iterated many times and then a fine search for the period will be found. The most important step in this case will be the function `delta_finder`. It chooses the best period from the outpus which are the PC anv V. The PC is transformed to a scalar version of it which is in the direction of the hyperdiagonal. Then the selection corresponds to the maximum scalar minus the average from the left scalars times the corresponding eigenvalue to the maximum chosen. Please refer to the code or the documentation to clarify this point.

### Outputs 

The principal outputs of the program will be the second and first periods from the iterations. It is possible to obtain other values from the residual data.

1. `period_final1`, `period_final2`, correspond to the best period found using the method. In an ascii file.
2. Eigenvalues, in a `.npy` output format.
3. Scalars (hyperdiagonal unitary vector times the eigenvector), in a `.npy` output format.
4. MSTEV (Maximum Scalar Times EigenValue), which is the selected merit function for the analysis. Same output as before, `.npy`.
