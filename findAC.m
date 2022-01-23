function [A,C,rcnd] = findAC(s,n,l,R,meth,tol,printw)
%FINDAC  Finds the system matrices A and C of a discrete-time system, given the
%        system order and the relevant part of the R factor of the concatenated
%        block-Hankel matrices, using subspace identification techniques (MOESP 
%        or N4SID).
% 
%        [A,C] = FINDAC(S,N,L,R,METH,TOL,PRINTW)  computes the system matrices
%        A and C. The model structure is:
%
%             x(k+1) = Ax(k) + Bu(k) + Ke(k),   k >= 1,
%             y(k)   = Cx(k) + Du(k) + e(k),
%
%        where x(k) and y(k) are vectors of length N and L, respectively.
%
%        [A,C,RCND] = FINDAC(S,N,L,R,METH,TOL,PRINTW)  also returns the vector
%        RCND of length 4 containing the condition numbers of the matrices
%        involved in rank decisions.
%
%        S is the number of block rows in the block-Hankel matrices.
%
%        METH is an option for the method to use:
%        METH = 1 :  MOESP method with past inputs and outputs;
%             = 2 :  N4SID method.
%        Default:    METH = 1.
%        Matrix R, computed by FINDR, should be determined with suitable arguments
%        METH and JOBD.
%
%        TOL is the tolerance used for estimating the rank of matrices. 
%        If  TOL > 0,  then the given value of  TOL  is used as a lower bound
%        for the reciprocal condition number.
%        Default:    prod(size(matrix))*epsilon_machine where epsilon_machine
%                    is the relative machine precision.
%
%        PRINTW is a switch for printing the warning messages.
%        PRINTW = 1: print warning messages;
%               = 0: do not print warning messages.
%        Default:    PRINTW = 0.
%
%        See also FINDABCD, FINDBD, FINDBDK, FINDR, ORDER, SIDENT
%

%        RELEASE 2.0 of SLICOT System Identification Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima 18-01-2000.
%
%        Revisions:
%        V. Sima 02-03-2009.
%   

nin = nargin;  nout = nargout;
% 
if nin < 7;  printw = 0;  end;
if nin < 6 || isempty(tol);    tol  = 0;  end;  
if nin < 5 || isempty(meth);   meth = 1;  end;  
if nin < 4, 
   error('Wrong number of input arguments')
end   
%
% Compute system matrices A and C.
job = 2;  nsmpl = 0;
if nout == 1,
   A = sident(meth,job,s,n,l,R,tol,nsmpl,[],[],printw);
elseif nout == 2,
   [A,C] = sident(meth,job,s,n,l,R,tol,nsmpl,[],[],printw);
elseif nout == 3,
   [A,C,rcnd] = sident(meth,job,s,n,l,R,tol,nsmpl,[],[],printw);
else
   error('Wrong number of output arguments')
end   
%
% end findAC
