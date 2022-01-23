function [B,C,D,x0] = o2s(n,m,l,theta,apply)
%O2S     Transforms a linear discrete-time system given in the
%        output normal form to a state-space representation.
%
%        [SYS,X0] = O2S(N,M,L,THETA,APPLY)  computes a state-space
%        representation of order N, SYS = (A,B,C,D) (an ss object),
%        with M inputs and L outputs, L > 0, with the initial state X0, 
%        given its output normal form with parameter vector THETA of 
%        length N*(L+M+1)+L*M.
%
%        APPLY is an optional integer specifying whether or not the
%        parameter vector THETA should be transformed using a 
%        bijective mapping:
%        APPLY = 1: apply the bijective mapping to the N vectors in
%                   THETA corresponding to the matrices A and C;
%        APPLY = 0: do not apply the bijective mapping.
%        Default:  APPLY = 0.
%        The value of APPLY should be the same as that used in the
%        corresponding call of the paired S2O function.
%
%        If L = 0, the matrix A cannot be recovered, and the function
%        should be called as follows:
%        [B,C,D,X0] = O2S(N,M,L,THETA,APPLY)
%
%        See also S2O
%

%        RELEASE 2.0 of SLICOT System Identification Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, Feb. 2002.
%
%        Revisions: Mar. 2009.
%

nin = nargin;
% 
if nin < 5,  
    apply = 0;
else
    if ~( numel( apply ) == 1 ),
        error( 'APPLY must be an integer variable' )
    end
    if ~( apply == 0 || apply == 1 ),
        apply = 0;  
    end
end;
%
if nin < 4, 
   error( 'Wrong number of input arguments' )
end   
%
if n < 0, 
   error( 'The system order should be non-negative' )
elseif m < 0, 
   error( 'The number of inputs should be non-negative' )
elseif l < 0, 
   error( 'The number of outputs should be non-negative' )
end   
%
if length( theta ) < n*( l + m + 1 ) + l*m, 
   error( 'The vector THETA is too short' )
end   
%
if l == 0,
    % System with no outputs. The matrix A cannot be recovered.
    [ B, C, D, x0 ] = onf2ss( n, m, l, theta, apply );
else
    % Below, C means x0 !
    [ As, Bs, Cs, Ds, C ] = onf2ss( n, m, l, theta, apply );
    % Below, B means sys !
    B = ss( As, Bs, Cs, Ds, 1);
end
%
% end o2s
