%SIDENT  MEX-function for computing a discrete-time state-space realization
%        (A,B,C,D) and Kalman gain K using SLICOT routine IB01BD.
%
%   [(A,C)(,B(,D))(,K,Q,Ry,S)(,rcnd)] = sident(meth,job,s,n,l,R(,tol,t,A,
%                                              C,printw))
%
%                 [A,C,B,D] = sident(meth,1,s,n,l,R)
%   [A,C,B,D,K,Q,Ry,S,rcnd] = sident(meth,1,s,n,l,R,tol,t)
%                     [A,C] = sident(meth,2,s,n,l,R)
%                         B = sident(meth,3,s,n,l,R,tol,0,A,C)
%         [B,K,Q,Ry,S,rcnd] = sident(meth,3,s,n,l,R,tol,t,A,C)
%                     [B,D] = sident(meth,4,s,n,l,R,tol,0,A,C)
%       [B,D,K,Q,Ry,S,rcnd] = sident(meth,4,s,n,l,R,tol,t,A,C)
%
%   SIDENT computes a state-space realization (A,B,C,D) and the Kalman
%   predictor gain K of a discrete-time system, given the system
%   order and the relevant part of the R factor of the concatenated 
%   block-Hankel matrices, using subspace identification techniques 
%   (MOESP, N4SID, or their combination).
%
%   The model structure is :
%
%         x(k+1) = Ax(k) + Bu(k) + Ke(k),   k >= 1,
%         y(k)   = Cx(k) + Du(k) + e(k),
%
%   where  x(k)  is the  n-dimensional state vector (at time k),
%          u(k)  is the  m-dimensional input vector,
%          y(k)  is the  l-dimensional output vector,
%          e(k)  is the  l-dimensional disturbance vector,
%   and  A, B, C, D, and K  are real matrices of appropriate dimensions.
%
%   Description of other input parameters:
%   meth  - integer option to determine the method to use:
%           = 1 : MOESP method with past inputs and outputs;
%           = 2 : N4SID method;
%           = 3 : combined method: A and C via MOESP, B and D via N4SID.
%   job   - integer option to determine the calculation to be performed:
%           = 1 : compute all system matrices, A, B, C, D;
%           = 2 : compute the matrices A and C only;
%           = 3 : compute the matrix B only;
%           = 4 : compute the matrices B and D only.
%   s     - the number of block rows in the processed input and output
%           block Hankel matrices.  s > 0.
%   R     - the 2*(m+l)*s-by-2*(m+l)*s part of  R  contains the
%           processed upper triangular factor  R  from the
%           QR factorization of the concatenated block-Hankel matrices,
%           and further details needed for computing system matrices.
%   tol   - (optional) tolerance used for estimating the rank of
%           matrices. If  tol > 0,  then the given value of  tol  is
%           used as a lower bound for the reciprocal condition number;
%           an m-by-n matrix whose estimated condition number is less
%           than  1/tol  is considered to be of full rank.
%           Default:    m*n*epsilon_machine where epsilon_machine is
%           the relative machine precision.
%   t     - (optional) the total number of samples used for calculating
%           the covariance matrices.  Either t = 0, or t >= 2*(m+l)*s.
%           This parameter is not needed if the covariance matrices
%           and/or the Kalman predictor gain matrix are not desired.
%           If t = 0, then K, Q, Ry, and S are not computed.
%           Default:    t = 0.
%   printw- (optional) switch for printing the warning messages.
%           = 1:  print warning messages;
%           = 0:  do not print warning messages.
%           Default:    printw = 0.
%
%   Description of other output parameters:
%   Q     - (optional) the n-by-n positive semidefinite state covariance
%           matrix used as state weighting matrix when computing the
%           Kalman gain.
%   Ry    - (optional) the l-by-l positive (semi)definite output
%           covariance matrix used as output weighting matrix when
%           computing the Kalman gain.
%   S     - (optional) the n-by-l state-output cross-covariance matrix
%           used as cross-weighting matrix when computing the Kalman
%           gain.
%   rcnd  - (optional) vector of length lr, containing estimates of the
%           reciprocal condition numbers of the matrices involved in
%           rank decisions, least squares, or Riccati equation solutions,
%           where lr = 4,  if Kalman gain matrix K is not required, and
%                 lr = 12, if Kalman gain matrix K is required.
%
%   Comments
%   1. The n-by-n system state matrix A, and the p-by-n system output 
%      matrix C are computed for job <= 2.
%   2. The n-by-m system input matrix B is computed for job <> 2.
%   3. The l-by-m system matrix D is computed for job = 1 or 4.
%   4. The n-by-l Kalman predictor gain matrix K and the covariance
%      matrices Q, Ry, and S are computed for t > 0.
%
% See also FINDBD, ORDER
%

% RELEASE 2.0 of SLICOT System Identification Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, Oct. 1999.
%
% Revisions:
%   V. Sima, May 2000, July 2000.
%
