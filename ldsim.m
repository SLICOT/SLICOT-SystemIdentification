% LDSIM   MEX-function for computing the output response of a 
%         linear discrete-time system using SLICOT routine TF01MY.
%
%   [Y(,x)] = ldsim(A,B,C,D,U(,x,ldwork))
%
%   LDSIM computes the output vector sequence y(1), y(2),..., y(t)
%
%        x(k+1) = A x(k) + B u(k)
%        y(k)   = C x(k) + D u(k),
%
%   given an initial state vector x(1), and the input vector sequence 
%   u(1), u(2),..., u(t), where y(k) and u(k) are vectors of length p
%   and m, respectively. The input trajectories are given as
%
%            ( u(1)' )
%            ( u(2)' )
%        U = (   :   ),  where ' denotes the transposition,
%            (   :   )
%            ( u(t)' )
%
%   and the output trajectories result in similarly.
%
%   Description of other input parameters:
%   x      - (optional) the initial state x(1).
%            Default: x = 0.
%   ldwork - (optional) the length of working array.
%            Default: ldwork = n.
%            Larger values could increase the efficiency.
%
%   Description of other output parameters:
%   x      - (optional) the final state x(t+1).
%
%   Comments
%   1. This function is faster than Matlab function lsim.
%
%   See also DSIM
%

%   RELEASE 2.0 of SLICOT System Identification Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
