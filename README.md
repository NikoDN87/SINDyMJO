# SINDyMJO
Implementation of Sparse Identification fo Nonlinear Dynamics Method for Madden-Julian Oscillation.

SINDyMJO.R is and R script code in which the method Sparse Identification fo Nonlinear Dynamics is implemented in order to infer two dimensional models of the Madden-Julian Oscillation.
The code uses OMI or RMM principal components as input data -- other simmilar indexes can be easly introduced.
The main idea of the code is to:
  1. Select random data of the MJO
  2. Apply SINDy (sequentially thresholded least square) to model the selected data
  3. Select models according to Akaike Information Criterion

The script presents the basic code used in: "Data Driven Models of the Madden-Julian Oscillation: understanding its evolution and ENSO modulation" by N.Diaz, N.Rubido and M.Barreiro.
