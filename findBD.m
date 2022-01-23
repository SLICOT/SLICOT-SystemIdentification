%FINDBD  MEX-function for estimating the initial state and the system
%        matrices B and D of a discrete-time system, using SLICOT routine 
%        IB01CD.
%
%   [(x0)(,B(,D))(,V)(,rcnd)] = findBD(jobx0,comuse(,job),A(,B),C(,D),Y
%                                      (,U,tol,printw,ldwork))
%
%     [x0,B,V,rcnd] = findBD(1,1,1,A,C,Y,U)
%   [x0,B,D,V,rcnd] = findBD(1,1,2,A,C,Y,U)
%        [B,V,rcnd] = findBD(2,1,1,A,C,Y,U)
%      [B,D,V,rcnd] = findBD(2,1,2,A,C,Y,U)
%       [x0,V,rcnd] = findBD(1,2,1,A,B,C,Y,U)
%       [x0,V,rcnd] = findBD(1,2,2,A,B,C,D,Y,U)
%         [x0,rcnd] = findBD(2,2)                % (Set x0 = 0, rcnd = 1)
%       [x0,V,rcnd] = findBD(1,3,A,C,Y)
%
%   Note: the example lines above may contain at the end the parameters
%         tol, printw, ldwork.
%
%   FINDBD estimates the initial state and/or the system matrices B and D
%   of a discrete-time system, given the system matrices A, C, and possibly
%   B, D, and the input and output trajectories of the system.
%
%   The model structure is :
%
%         x(k+1) = Ax(k) + Bu(k),   k >= 1,
%         y(k)   = Cx(k) + Du(k),
%
%   where  x(k)  is the  n-dimensional state vector (at time k),
%          u(k)  is the  m-dimensional input vector,
%          y(k)  is the  l-dimensional output vector,
%   and  A, B, C, and D  are real matrices of appropriate dimensions.
%
%   Description of other input parameters:
%   jobx0 - integer option to specify whether or not the initial state 
%           should be computed:
%           = 1 : compute the initial state x(0);
%           = 2 : do not compute the initial state (possibly, because
%                 x(0) is known to be zero).
%   comuse- integer option to specify whether the system matrices B
%           and D should be computed or used:
%           = 1 : compute the matrices B and D, as specified by job;
%           = 2 : use the matrices B and D, as specified by job;
%           = 3 : do not compute/use the matrices B and D.
%   job   - integer option to determine which of the system matrices
%           B and D should be computed or used:
%           = 1 : compute/use the matrix B only (D is known to be zero);
%           = 2 : compute/use the matrices B and D.
%           job must not be specified if jobx0 = 2 and comuse = 2, or
%           if comuse = 3.
%   Y     - the t-by-l output-data sequence matrix.  Column  j  of  Y 
%           contains the  t  values of the j-th output component for
%           consecutive time increments.
%   U     - the t-by-m input-data sequence matrix (input when jobx0 = 1 
%           and comuse = 2, or comuse = 1).  Column  j  of  U 
%           contains the  t  values of the j-th input component for
%           consecutive time increments.
%   tol   - (optional) tolerance used for estimating the rank of
%           matrices. If  tol > 0,  then the given value of  tol  is
%           used as a lower bound for the reciprocal condition number;
%           an m-by-n matrix whose estimated condition number is less
%           than  1/tol  is considered to be of full rank.
%           Default:    m*n*epsilon_machine where epsilon_machine is the
%           relative machine precision.
%   printw- (optional) switch for printing the warning messages.
%           = 1:  print warning messages;
%           = 0:  do not print warning messages.
%           Default:    printw = 0.
%   ldwork- (optional) the workspace size.
%           Default :   computed by the formula
%           LDWORK = MAX( minimum workspace size needed, 2*CSIZE/3,
%                         CSIZE - ( m + l )*t - 2*n*( n + m + l ) - l*m )
%           where CSIZE is the cache size in double precision words.
%
%   Description of other output parameters:
%   V     - the n-by-n orthogonal matrix which reduces A to a real Schur
%           form (output when jobx0 = 1 or comuse = 1).
%   rcnd  - (optional) the reciprocal condition numbers of the matrices
%           involved in rank decisions.
%
%   Comments
%   1. The n-by-m system input matrix B is an input parameter when jobx0 = 1 
%      and comuse = 2, and it is an output parameter when comuse = 1.
%   2. The l-by-m system matrix D is an input parameter when jobx0 = 1, 
%      comuse = 2 and job = 2, and it is an output parameter when comuse = 1 
%      and job = 2.
%   3. The n-vector of estimated initial state x(0) is an output parameter
%      when jobx0 = 1, but also when jobx0 = 2 and comuse <= 2, in which case
%      it is set to 0.
%   4. If ldwork is specified, but it is less than the minimum workspace size 
%      needed, that minimum value is used instead.
%
% See also INISTATE, FINDX0BD
%

%   RELEASE 2.0 of SLICOT System Identification Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   Contributor:
%   V. Sima, Katholieke Univ. Leuven, Belgium, May 2000.
%
%   Revisions:
%   V. Sima, July 2000.
%
