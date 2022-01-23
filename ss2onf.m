%SS2ONF  MEX-function for converting a linear discrete-time system
%        into the output normal form using SLICOT routine TB01VD.
%
%   theta = ss2onf(A,B,C,D,x0,apply)
%
%   SS2ONF transforms a linear discrete-time system given as (A,B,C,D),
%   with initial state x0, into the output normal form, with parameter
%   vector theta. The matrix A is assumed to be stable.
%   Optionally, a bijective transformation is used, guaranteeing that
%   the inverse conversion (using ONF2SS) can always be performed.
%
%   Description of other input parameters:
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
%
% See also ONF2SS
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
