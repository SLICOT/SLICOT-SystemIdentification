# **SLICOT System Identification Toolbox**  

## About 

The `SLICOT System Identification Toolbox ` (`SLICOT-SystemIdentification`) includes [SLICOT](http://slicot.org/)-based MATLAB and Fortran tools for linear and Wiener-type, time-invariant discrete-time multivariable systems. Subspace-based approaches MOESP - Multivariable Output-Error state SPace identification, N4SID - Numerical algorithms for Subspace State Space System IDentification, and their combination, are used to identify linear systems, and to initialize the parameters of the linear part of a Wiener system. All parameters of a Wiener system are then estimated using a specialized Levenberg-Marquardt algorithm.

The main functionalities of the toolbox include:

    * identification of linear discrete-time state space systems (A, B, C, D)
    * identification of state and output (cross-)covariance matrices for such systems
    * estimation of the associated Kalman gain matrix
    * estimation of the initial state
    * conversion from/to a state-space representation to/from the output normal form parameterization
    * identification of discrete-time Wiener systems
    * computation of the output response of Wiener systems.

The toolbox main features are:

    * computational reliability
    * high numerical efficiency, using structure exploiting algorithms and dedicated linear algebra tools
    * possible speed-up factors larger then 10 in comparison with the commonly used software tools flexibility and easy-of-use
    * ability to process multiple (possibly connected) data batches
    * standardized interfaces

The programs have been extensively tested on various test examples and are fully documented.

The current release of `SLICOT-SystemIdentification` is version 1.0, dated November 1, 2021.

## Requirements

The codes have been tested with MATLAB 2015b through 2021b. To use the functions, the Control System Toolbox must be installed in MATLAB running under 64-bit Windows 7, 8, 8.1 or 10.  

## License

* See [`LICENSE`](https://github.com/SLICOT/SLICOT-BasicControl/blob/master/LICENSE) for licensing information.

## References

Please cite `SLICOT-SystemIdentification` using at least one of the following references: 

* P. Benner, D. Kressner, V, Sima, and A. Varga, [The SLICOT Toolboxes - a Survey](http://slicot.org/objects/software/reports/SLWN2009-1.pdf), _SLICOT Working Note 2009-1, August 2009._
* P. Benner, D. Kressner, V. Sima, A. Varga, Die SLICOT-Toolboxen für Matlab - The SLICOT Toolboxes for Matlab (in German), _at – Automatisierungstechnik, 58 (2010)._

