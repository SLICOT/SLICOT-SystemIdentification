% ORDER.F  - MEX-function for computing the order of a discrete-time
%            system using SLICOT routine IB01AD.
%
%   [R(,n,sval,rcnd)] = order(meth,alg,jobd,batch,conct,s,Y(,U,tol,
%                             printw,ldwork,R))
%
%   For one block (data sequences Y, U):
%     [R,n,sval,rcnd] = order(meth,alg,jobd,4,conct,s,Y,U); 
%
%   For f blocks (data sequences Yj, Uj, j = 1 : f):
%   R = order(meth,alg,jobd,1,conct,s,Y1,U1);
%   for j = 2 : f - 1
%      R = order(meth,alg,jobd,2,conct,s,Yj,Uj,tol,printw,ldwork,R)
%   end
%   [R,n,sval,rcnd] = order(meth,alg,jobd,3,conct,s,Yf,Uf,tol);
%
%   ORDER preprocesses the input-output data for estimating the matrices 
%   of a linear time-invariant dynamical system, using Cholesky or (fast)
%   QR factorization and subspace identification techniques (MOESP and
%   N4SID), and then estimates the order of a discrete-time realization.
%
%   The model structure is :
%
%         x(k+1) = Ax(k) + Bu(k) + w(k),   k >= 1,
%         y(k)   = Cx(k) + Du(k) + e(k),
%
%   where  x(k)  is the  n-dimensional state vector (at time k),
%          u(k)  is the  m-dimensional input vector,
%          y(k)  is the  l-dimensional output vector,
%          w(k)  is the  n-dimensional state disturbance vector,
%          e(k)  is the  l-dimensional output disturbance vector,
%   and  A, B, C, and D  are real matrices of appropriate dimensions.
%
%   Description of input parameters:
%   meth  - integer option to determine the method to use:
%           = 1 : MOESP method with past inputs and outputs;
%           = 2 : N4SID method.
%   alg   - integer option to determine the algorithm for computing the
%           triangular factor of the concatenated block-Hankel matrices
%           built from the input-output data:
%           = 1 : Cholesky algorithm on the correlation matrix;
%           = 2 : fast QR algorithm;
%           = 3 : standard QR algorithm.
%   jobd  - integer option to specify if the matrices B and D should
%           later be computed using the MOESP approach:
%           = 1 : the matrices B and D should later be computed using
%                 the MOESP approach;
%           = 2 : the matrices B and D should not be computed using
%                 the MOESP approach.
%           This parameter is not relevant for meth = 2.
%   batch - integer option to specify whether or not sequential data
%           processing is to be used, and, for sequential processing,
%           whether or not the current data block is the first block,
%           an intermediate block, or the last block, as follows:
%           = 1 : the first block in sequential data processing;
%           = 2 : an intermediate block in sequential data processing;
%           = 3 : the last block in sequential data processing;
%           = 4 : one block only (non-sequential data processing).
%   conct - integer option to specify whether or not the successive data
%           blocks in sequential data processing belong to a single
%           experiment, as follows:
%           = 1 : the current data block is a continuation of the
%                 previous data block and/or it will be continued by the
%                 next data block;
%           = 2 : there is no connection between the current data block
%                 and the previous and/or the next ones.
%           This parameter is not used if batch = 4.
%   s     - the number of block rows in the input and output block
%           Hankel matrices to be processed.  s > 0.
%   Y     - the t-by-l output-data sequence matrix.  Column j of  Y
%           contains the  t  values of the j-th output component for
%           consecutive time increments.
%   U     - (optional) the t-by-m input-data sequence matrix.  Column j
%           of  U  contains the  t  values of the j-th input component
%           for consecutive time increments.
%           Default:    U = [].
%   tol   - (optional) vector of length 2 containing tolerances:
%           tol(1) - tolerance used for estimating the rank of matrices.
%           If  tol(1) > 0,  then the given value of  tol(1)  is used
%           as a lower bound for the reciprocal condition number;
%           an m-by-n matrix whose estimated condition number is less
%           than  1/tol(1)  is considered to be of full rank.
%           If  tol(1) <= 0,  then a default value m*n*epsilon_machine
%           is used, where epsilon_machine is the relative machine
%           precision.
%           tol(2) - tolerance used for determining an estimate of the
%           system order. If  tol(2) >= 0,  the estimate is indicated
%           by the index of the last singular value greater than or
%           equal to  tol(2).  (Singular values less than  tol(2)  are
%           considered as zero.) When  tol(2) = 0,  an internally
%           computed default value,  tol(2) = s*epsilon_machine*sval(1),
%           is used, where  sval(1)  is the maximal singular value, and
%           epsilon_machine the relative machine precision.
%           When  tol(2) < 0,  the estimate is indicated by the index of
%           the singular value that has the largest logarithmic gap to
%           its successor.
%           Default:    tol(1:2) = [0,-1].
%   printw- (optional) switch for printing the warning messages.
%           = 1:  print warning messages;
%           = 0:  do not print warning messages.
%           Default:    printw = 0.
%   ldwork- (optional) the workspace size.
%           Default : computed by the formulas
%           nr = 2*( m + l )*s
%           LDWORK = ( t - 2*s + 3 + 64 )*nr
%           if ( CSIZE > MAX( nr*nr + t*( m + l ) + 16, 2*nr ) ) then
%              LDWORK = MIN( LDWORK, CSIZE - nr*nr - t*( m + l ) - 16 )
%           else
%              LDWORK = MIN( LDWORK, MAX( 2*nr, CSIZE/2 ) )
%           end if
%           LDWORK = MAX( minimum workspace size needed, LDWORK )
%           where CSIZE is the cache size in double precision words.
%           If LDWORK is specified less than the minimum workspace size 
%           needed, that minimum value is used instead.
%   R     - (optional) if batch = 2 or 3, the 2*(m+l)*s-by-2*(m+l)*s
%           (upper triangular, if alg <> 2) part of  R  must contain the
%           (upper triangular) matrix  R  computed at the previous call
%           of this mexfile in sequential data processing. If conct = 1,
%           R  has an additional column, also set at the previous call. 
%           If alg = 2,  R  has m+l+1 additional columns, set at the
%           previous call. 
%           This parameter is not used for batch = 1 or batch = 4.
%
%   Description of output parameters:
%   R     - if batch = 3 or 4, the 2*(m+l)*s-by-2*(m+l)*s part of  R 
%           contains the processed upper triangular factor  R  from the
%           QR factorization of the concatenated block-Hankel matrices,
%           and further details needed for computing system matrices.
%           If batch = 1 or 2, then  R  contains intermediate results
%           needed at the next call of this mexfile. If batch = 1 or 2
%           and conct = 1,  R  has an additional column, also set before
%           return. If batch = 1 or 2 and alg = 2,  R  has m+l+1 
%           additional columns, set before return.
%   n     - the order of the system.
%   sval  - (optional) the singular values used for estimating the order
%           of the system.
%   rcnd  - (optional) if meth = 2, vector of length 2 containing the
%           reciprocal condition numbers of the matrices involved in
%           rank decisions or least squares solutions.
%
%   Comments
%   1. The Cholesy or fast QR algorithms can be much faster (for large data
%      blocks) than QR algorithm, but they cannot be used if the correlation 
%      matrix, H'*H, is not positive definite. In such a case, the code
%      automatically switches to the QR algorithm, if sufficient workspace
%      is provided and batch = 4. 
%   2. If ldwork is specified, but it is less than the minimum workspace size 
%      needed, that minimum value is used instead.
%
% See also FINDBD, SIDENT
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
