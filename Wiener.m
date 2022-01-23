%WIENER  MEX-function for computing the output of a Wiener 
%        system using SLICOT routine NF01AD.
%
%   y = Wiener(n,l,nn,theta,u(,ldwork))
%
%   WIENER computes the output y of the Wiener system
%
%        x(k+1) = A*x(k) + B*u(k)
%        z(k)   = C*x(k) + D*u(k),
%
%        y(k)   = f(z(k),wb(1:l)),
%
%   where the linear discrete-time system is given as its output normal
%   form, with parameter vector theta(l1+1:l1+l2), and the parameters of
%   the nonlinear part are contained in theta(1:l1), with
%   l1 = (nn*(l+2)+1)*l, l2 = n*(l+m+1)+l*m.
%
%   Description of input parameters:
%   n      - the order of the linear system.
%   l      - the number of the system outputs.
%   nn     - the number of neurons of the nonlinear part.
%   theta  - the (nn*(l+2)+1)*l+n*(l+m+1)+l*m parameter vector.
%   u      - the t-by-m input trajectory.
%   ldwork - (optional) the length of working array.
%            Default: ldwork = t*l+ max( 2*nn,(n+l)*(n+m)+2*n+
%                                             max(n*(n+l),nml)),
%                              where nml = n+m+l, if m > 0, and
%                                    nml = l,     if m = 0.
%            Larger values could increase the efficiency.
%
%   Description of output parameters:
%   y      - the t-by-l output trajectory.
%
% See also wident
%

% RELEASE 2.0 of SLICOT System Identification Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, Apr. 2001.
%
% Revisions:
%   Feb. 2002.
%
