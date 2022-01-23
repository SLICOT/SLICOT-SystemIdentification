%ONF2SS  MEX-function for converting a linear discrete-time system
%        given in the output normal form to a state-space
%        representation using SLICOT routine TB01VY.
%
%   [A,B,C,D,x0] = onf2ss(n,m,l,theta,apply)  % if l > 0;
%     [B,C,D,x0] = onf2ss(n,m,l,theta,apply)  % if l = 0.
%
%   ONF2SS transforms a linear discrete-time system given as its output
%   normal form, with parameter vector theta, into the state-space
%   representation (A,B,C,D), with the initial state x0.
%
%   Description of other input parameters:
%   n      - the order of the system.
%   m      - the number of the system inputs.
%   l      - the number of the system outputs.
%   apply  - (optional) integer specifying whether or not the parameter
%            vector should be transformed using a bijective mapping:
%            = 1: apply the bijective mapping to the n vectors in
%                 theta corresponding to the matrices A and C;
%            = 0: do not apply the bijective mapping.
%            Default:  apply = 0.
%
%   Comments
%   1. The vector theta has n*(l+m+1)+l*m elements, where [l,m] = size(D),
%      and n is the order of the matrix A.
%   2. The value of apply should be the same as that used in the
%      corresponding call of the paired SS2ONF function.
%
% See also SS2ONF
%

% RELEASE 2.0 of SLICOT System Identification Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributors:
%   A. Riedel, R. Schneider, Chemnitz Univ. of Technology, Mar. 2001.
%
% Revisions:
%   V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001,
%   Feb. 2002.
%
