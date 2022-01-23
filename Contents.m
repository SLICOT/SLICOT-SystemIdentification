% SLIDENT - SLICOT System Identification Toolbox.
% RELEASE 2.0 of SLICOT Toolboxes.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% This toolbox includes M- and MEX-functions for system identification of
% linear time-invariant state-space systems and nonlinear Wiener systems.
%
% MATLAB M-files (except help files).
% -----------------------------------
%
%  Linear time-invariant multivariable state-space systems identification
%  using subspace techniques.
%
%   a. Method-oriented functions (preferred choice: slmoen4):
%   slmoen4     - Finds the system matrices and the Kalman gain of a discrete-time
%                 system, using combined MOESP and N4SID subspace identification
%                 techniques.
%   slmoesm     - Finds the system matrices, the Kalman gain, and initial state of a
%                 discrete-time system, using combined MOESP subspace identification
%                 and system simulation techniques.
%   slmoesp     - Finds the system matrices and the Kalman gain of a discrete-time
%                 system, using MOESP subspace identification technique.
%   sln4sid     - Finds the system matrices and the Kalman gain of a discrete-time
%                 system, using N4SID subspace identification technique.
%
%   b. Lower-level functions:
%   findR       - Preprocesses the input-output data for estimating the matrices
%                 of a linear time-invariant dynamical system, using Cholesky or
%                 (fast) QR factorization and subspace identification techniques 
%                 (MOESP or N4SID), and estimates the system order.
%   findABCD    - Finds the system matrices and the Kalman gain of a discrete-time
%                 system, given the system order and the relevant part of the
%                 R factor of the concatenated block-Hankel matrices.
%   findAC      - Finds the system matrices A and C of a discrete-time system,
%                 given the system order and the relevant part of the R factor
%                 of the concatenated block-Hankel matrices.
%   findBDK     - Finds the system matrices B and D and the Kalman gain of a
%                 discrete-time system, given the system order, the matrices A
%                 and C, and the relevant part of the R factor of the
%                 concatenated block-Hankel matrices.
%   findx0BD    - Estimates the initial state and/or the matrices B and D of a
%                 discrete-time linear system, given the (estimated) system
%                 matrices A, C, and a set of input/output data.
%   inistate    - Estimates the initial state of a discrete-time system, given the
%                 (estimated) system matrices, and a set of input/output data.
%   dsim        - Computes the output response of a linear discrete-time system.
%
%  Nonlinear Wiener systems identification.
%
%   NNout       - Computes the output of a set of neural networks used to model
%                 the nonlinear part of a Wiener system.
%   o2s         - Transforms a linear discrete-time system given in the
%                 output normal form to a state-space representation.
%   s2o         - Transforms a state-space representation of a linear
%                 discrete-time system into the output normal form.
%
% MATLAB MEX-files.
% -----------------
%
%  Linear time-invariant multivariable state-space systems identification
%  using subspace techniques.
%
%   order       - Preprocesses the input-output data and estimate the
%                 system order.
%   sident      - Estimates the system matrices, the covariance matrices and
%                 Kalman gain.
%   findBD      - Estimates the initial condition and/or the matrices B and D
%                 using the (other) system matrices and the input-output data.
%
%  Nonlinear Wiener systems identification.
%
%   Wiener      - Computes the output of a Wiener system.
%   onf2ss      - Transforms a linear discrete-time system given in the
%                 output normal form to a state-space representation.
%   ss2onf      - Transforms a state-space representation of a linear
%                 discrete-time system into the output normal form.
%   wident      - Computes a discrete-time model of a Wiener system
%                 using a neural network approach and a MINPACK-like
%                 Levenberg-Marquardt algorithm.
%   widentc     - Computes a discrete-time model of a Wiener system
%                 using a neural network approach and a Levenberg-
%                 Marquardt algorithm (with a Cholesky-based, or
%                 a conjugate gradients solver).
%
%  Auxiliary MEX-file.
%
%   ldsim       - Computes the output response of a linear discrete-time 
%                 system (faster than the MATLAB function lsim).
%
%  Note: the MEX-files wident and widentc have to be called directly, 
%        since no M-file is available.
%
% Demonstrations.
% ---------------
%
%   slidemo     - Linear System Identification Demonstration.
%   slwidemo    - Wiener System Identification Demonstration.
%   cmp_ident   - Comparison between SLIDENT function slmoen4 and
%                 MATLAB System Identification function n4sid.
%   test_dsim   - Comparison between SLIDENT function dsim and MATLAB
%                 function lsim.
%
% Note:
% -----
%
%   Command-line help functions are available for all MATLAB M- and 
%   MEX-functions included in this toolbox.
%

%  ..CONTRIBUTOR..
%
%   V. Sima, Katholieke Univ. Leuven, Belgium, Sept. 2001.
%
%   Revisions:
%   V. Sima,  April 14, 2002, August 31, 2002, March 15, 2003,
%             May 15, 2003, August 25, 2003, April 11, 2004, April 20, 2004,
%             January 25, 2005, May 25, 2005, March 9, 2009, Jan. 2022.




