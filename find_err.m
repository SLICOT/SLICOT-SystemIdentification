function [err,ye,x0] = find_err(y,u,sys,K)
%FIND_ERR  Computes the estimated output and its relative error, for
%          a linear discrete-time system, using the input and output
%          trajectories, and the identified system matrices and
%          Kalman gain.
%
%          ERR = FIND_ERR(Y,U,SYS)  computes the relative error,
%          norm(Y - Ye,1)/norm(Y,1), where Ye is determined using
%          SYS = (A,B,C,D) (an ss object), the given input U, and
%          an estimate of the initial state, X0.
%
%          [ERR,Ye] = FIND_ERR(Y,U,SYS)  also returns Ye.
%          [ERR,Ye,X0] = FIND_ERR(Y,U,SYS)  also returns X0.
%
%          ERR = FIND_ERR(Y,U,SYS,K)  computes the relative error
%          above, using also the Kalman gain matrix K.
%
%          [ERR,Ye] = FIND_ERR(Y,U,SYS,K)  also returns Ye.
%          [ERR,Ye,X0] = FIND_ERR(Y,U,SYS,K)  also returns X0.

%        RELEASE 2.0 of SLICOT System Identification Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima 16-09-2000.
%
%        Revisions:
%        V. Sima 30-12-2000, 04-03-2009.
%

nin = nargin;
%
x0 = inistate(sys,y,u);
if nin == 4,
   x  = ltitr( ( sys.a - K*sys.c ),[ ( sys.b - K*sys.d ) K],[u y],x0);
else
   x  = ltitr(sys.a,sys.b,u,x0);
end
ye = x * sys.c.' + u * sys.d.';
%
err = norm(y - ye,1)/norm(y,1);
%
% end find_err
