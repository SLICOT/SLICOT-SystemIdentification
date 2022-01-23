%WIDENT  MEX-function for computing a discrete-time model of a
%        Wiener system using SLICOT routine IB03BD.
%
%   [xopt(,perf,nf,rcnd)] = wident(job,u,y,nn(,s,n,x,iter,nprint,tol,
%                                  seed,printw,ldwork))
%
%   [xopt(,perf,nf,rcnd)] = wident(1,u,y,nn,s,n,x(,iter,nprint,tol,
%                                  printw,ldwork))
%   [xopt(,perf,nf)]      = wident(2,u,y,nn,n,x(,iter,nprint,tol,
%                                  seed,printw,ldwork))
%   [xopt(,perf,nf,rcnd)] = wident(3,u,y,nn,s(,n,[],iter,nprint,tol,
%                                  seed,printw,ldwork))
%   [xopt(,perf,nf)]      = wident(4,u,y,nn,n,x(,iter,nprint,tol,
%                                  printw,ldwork))
%
%   WIDENT computes a set of parameters for approximating a Wiener system
%   in a least-squares sense, using a neural network approach and a
%   MINPACK-like Levenberg-Marquardt algorithm.
%
%   Description of input parameters:
%   job    - integer specifying which parameters must be initialized:
%            = 1 : initialize the linear part only;
%            = 2 : initialize the static nonlinearity only;
%            = 3 : initialize both linear and nonlinear parts;
%            = 4 : do not initialize anything, x already contains
%                  an initial approximation.
%   u      - the t-by-m input  trajectory, where t is the number of
%            samples.
%   y      - the t-by-l output trajectory.
%   nn     - the number of neurons used to approximate the nonlinear
%            part.
%   s      - if job = 1 or 3, the number of block rows in the input and
%            output block Hankel matrices to be processed.  s > 0.
%            The argument s must not be specified if job = 2 or 4.
%   n      - (optional) the order of the linear part of the system or
%            an option for how to compute it.
%            If job = 1 or 3, and n < 0, the order will be found by the
%            program. Otherwise, n should be non-negative (or positive,
%            but less than s, if job = 1 or 3).
%            The value of n should be given if job = 2 or 4.
%            Default for job = 1 or 3: the order found automatically.
%   x      - the lx0 vector of initial parameters, where
%            lx0 = l1,  if job = 1,
%            lx0 = l2,  if job = 2,
%            lx0 = lx,  if job = 4,
%            with l1 = (nn*(l+2)+1)*l, l2 = n*(l+m+1)+l*m, lx = l2+l1.
%            Here, l1 is the number of parameters for the nonlinear part
%            and l2 is the number of parameters for the linear part.
%            If n < 0 on entry, then s-1 is used instead of n for
%            setting l2 and lx, and for allocating the workspace.
%            If job = 3, the vector x is not needed.
%   iter   - (optional) vector of length 2 containing the maximal
%            numbers of iterations for the Levenberg-Marquardt algorithm
%            for the initialization of the static nonlinearity (ignored
%            if job = 1 or 4), and for the optimization process; these
%            numbers are stored in iter(1) and iter(2), respectively.
%            Default: iter = [10*lx 10*lx].
%   nprint - (optional) integer specifying the frequency of printing
%            details of the optimization process. If nprint > 0, the
%            intermediate results are printed in the file IB03BD.prn,
%            which is overwritten at each execution of this mexfile.
%            Default: nprint = 0 (no printing).
%   tol    - (optional) vector of length 2 containing the absolute
%            tolerances for the initialization of the static
%            nonlinearity (ignored if job = 1 or 4), and for the
%            optimization process, respectively.
%            Default: tol = [sqrt(eps) sqrt(eps)], with
%            eps = epsilon_machine.
%   seed   - (optional) if job = 2 or 3, vector of length 4 containing
%            the random number generator seed used to initialize the
%            parameters of the static nonlinearity.
%   printw - (optional) switch for printing the warning messages.
%            = 1:  print warning messages;
%            = 0:  do not print warning messages.
%            Default:    printw = 0.
%   ldwork - (optional) the length of working array.
%            Default: the minimum workspace needed, computed by the
%                     program based on the specified input arguments.
%            Larger values could increase the efficiency.
%
%   Description of output parameters:
%   xopt   - the lx vector of optimal parameters.
%            The corresponding linear system and its initial state can
%            be found as follows:
%            m = size(u,2);  l = size(y,2);  l1 = (nn*(l+2)+1)*l;
%            n = ( length(xopt) - l1 - l*m ) / ( l+m+1 );
%            [A,B,C,D,x0] = onf2ss(n,m,l,xopt(l1+1:end),1);
%   perf   - (optional) vector of length 4 or 8 containing performance
%            results:
%            perf(1) contains the optimal value of the workspace length;
%            perf(2) contains the maximum residual error norm;
%            perf(3) contains the total number of iterations performed;
%            perf(4) contains the final Levenberg factor.
%            If job = 2 or 3, then similar results for the nonlinear
%            part initialization step are returned:
%            perf(5) contains the optimal value of the workspace length;
%            perf(6) contains the maximum residual error norm;
%            perf(7) contains the total number of iterations performed;
%            perf(8) contains the final Levenberg factor (maximum over
%                    all outputs).
%   nf     - (optional) 2-vector containing additional performance
%            results:
%            nf(1) contains the (total) number of function evaluations;
%            nf(2) contains the (total) number of Jacobian evaluations.
%   rcnd   - (optional) if job = 1 or 3, vector of suitable length
%            containing the reciprocal condition number estimates for
%            determining the linear part by susbspace techniques.
%
% See also ONF2SS, SS2ONF
%

% RELEASE 2.0 of SLICOT System Identification Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, Apr. 2001.
%
% Revisions:
%   V. Sima, Katholieke University Leuven, Belgium, Feb. 2002.
%
