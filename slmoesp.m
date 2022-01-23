function [sys,K,rcnd,R] = slmoesp(s,y,u,n,R,tol,printw)
%SLMOESP  Finds the system matrices and the Kalman gain of a discrete-time 
%         system, using MOESP subspace identification technique.
% 
%        SYS = SLMOESP(S,Y,U,N,ALG,TOL,PRINTW)  computes a state-space
%        realization of order N, SYS = (A,B,C,D) (an ss object), from the
%        input and output data U and Y. The model structure is:
%
%             x(k+1) = Ax(k) + Bu(k) + Ke(k),   k >= 1,
%             y(k)   = Cx(k) + Du(k) + e(k),
%
%        where x(k) is a vector of length N, and the vectors y(k) and u(k) are
%        transposes of the k-th rows of Y and U matrices, respectively.
%        If N = 0, or N = [], or N is omitted from the input parameters, the user
%        is prompted to provide its value, after inspecting the singular values.
%        If N < 0, then N is determined automatically, according to TOL(2).
%
%        [SYS,K,RCND] = SLMOESP(S,Y,U,N,ALG,TOL,PRINTW)  also returns the Kalman
%        predictor gain, as well as the vector RCND of length 12 containing
%        estimates of accuracy parameters, mainly of the reciprocal condition
%        numbers of the matrices involved in rank decisions, least squares or 
%        Riccati equation solutions. RCND(12) is an estimate of the forward
%        error bound of the Riccati equation solution used for finding K.
%
%        S is the number of block rows in the block-Hankel matrices built by
%        the MOESP identification technique (S > N).
%
%        ALG is an option for the algorithm to compute the triangular factor of
%        the concatenated block-Hankel matrices built from the input-output data:
%        ALG = 1 :   Cholesky algorithm on the correlation matrix;
%            = 2 :   fast QR algorithm;
%            = 3 :   standard QR algorithm.
%        Default:    ALG = 1.
%
%        TOL is a vector of length 2 containing tolerances: 
%        TOL(1) is the tolerance for estimating the rank of matrices.
%        If  TOL(1) > 0,  the given value of  TOL(1)  is used as a
%        lower bound for the reciprocal condition number.
%        Default:    TOL(1) = prod(size(matrix))*epsilon_machine where
%                    epsilon_machine is the relative machine precision.
%        TOL(2) is the tolerance for estimating the system order.
%        If  TOL(2) >= 0,  the estimate is indicated by the index of
%        the last singular value greater than or equal to  TOL(2). 
%        (Singular values less than  TOL(2) are considered as zero.)
%        When  TOL(2) = 0,  then  S*epsilon_machine*sval(1)  is used instead
%        TOL(2),  where  sval(1)  is the maximal singular value.
%        When  TOL(2) < 0,  the estimate is indicated by the index of the
%        singular value that has the largest logarithmic gap to its successor.
%        Default:    TOL(2) = -1.
%
%        PRINTW is a switch for printing the warning messages.
%        PRINTW = 1: print warning messages;
%               = 0: do not print warning messages.
%        Default:    PRINTW = 0.
%
%        [SYS,K] = SLMOESP(S,Y)  identifies a stochastic system without inputs,
%        the order N being prompted to be set by the user.
%        [SYS,K] = SLMOESP(S,Y,[],N)  identifies a stochastic system, the 
%        order N being specified.
%
%        [SYS,K,RCND,R] = SLMOESP(S,Y,U,N,ALG,TOL,PRINTW)  also returns the 
%        processed upper triangular factor  R  of the concatenated block-Hankel 
%        matrices built from the input-output data. It can be used for fast
%        identification of systems of various orders, using, for instance, the
%        following commands:
%
%        [SYS,K,RCND,R] = SLMOESP(S,Y,U,N0,ALG,TOL,PRINTW);
%        for N = N0+1 : min( N0+Nf, S-1 )
%          [SYS,K,RCND] = SLMOESP(S,Y,U,N,R,TOL,PRINTW);
%        end
%        The data values for Y and U are not used inside the loop (only the
%        size of Y is needed), but R replaces ALG. Clearly, the systems of
%        various orders (from N0+1 to min( N0+NF, S-1 )), should be used 
%        inside the loop. 
%
%        See also ORDER, SIDENT, SLMOEN4, SLMOESM, SLN4SID
%

%        RELEASE 2.0 of SLICOT System Identification Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima 06-07-2000.
%
%        Revisions:
%        V. Sima, Nov. 8, 2004, Mar. 3, 2009.
%   

nin = nargin;  nout = nargout;
% 
if nin < 7,  printw = 0;  end;
if nin < 6 || isempty(tol),  tol(1:2) = [0,-1];  end;  
if nin < 5 || isempty(R),    R = 1;  end;  
if nin < 4 || isempty(n),    n = 0;  end; 
[t,l] = size(y); 
if nin < 3 || isempty(u), 
   u = zeros(t,0);
end   
job  = 1;   
jobd = job;
if nin < 2, 
   error('Wrong number of input arguments')
end   
if n >= s, 
   error('The system order N must be less than the number of block rows S')
end   
% Use MOESP.
meth = 1; 
%
nsmpl = t - 2*s + 1;
if nout == 1,  nsmpl = 0;  end
%
% Compute R and the system order.
if numel( R ) == 1,
   % Find R and the system order N. The input parameter R is here alg.
   [R,na,sval] = findR(s,y,u,meth,R,jobd,tol,printw);
   if n < 0,
      % N is determined automatically.
      n = min( na, s - 1);
   elseif n == 0,
      % N is set by the user.
      figure(gcf);  hold off;  subplot;
      h = bar( 1 : l*s, sval );
      a = get( h );  xx = a.XData;  yy = a.YData;
      semilogy( xx, yy + 10^( floor( log10( min( sval ) ) ) ) );
      axis( [0, length( sval ) + 1, 10^( floor( log10( min( sval ) ) ) ),...
                                    10^( ceil(  log10( max( sval) ) ) ) ] );
      title('Singular Values');
      xlabel('Order');
      disp(' ');
      disp(['System order defined by the chosen tolerance, N = ',sprintf('  %d',na),'.']);
      answer = input('Do you want to change it? Type y or yes (default), or n or no: ','s');
      if isempty(answer) || strcmp(answer,'y')   || strcmp(answer,'Y') || ...
                            strcmp(answer,'yes') || strcmp(answer,'YES'),
         disp(' ');
         n = 0;
         while ( n <= 0 ) || ( n >= s )
            n = input('System order ? ');
            if ( isempty(n) );  n = 0;  end
         end
      else
         n = min( na, s - 1);
      end
   end
end
%      
% Compute all system matrices.
if nout == 1,
   [A,C,B,D] = sident(meth,job,s,n,l,R,tol(1),nsmpl,[],[],printw);
elseif nout == 2,
   [A,C,B,D,K] = sident(meth,job,s,n,l,R,tol(1),nsmpl,[],[],printw);
elseif nout == 3 || nout == 4,
   [A,C,B,D,K,Q,Ry,S,rcnd] = sident(meth,job,s,n,l,R,tol(1),nsmpl,[],[],printw);
elseif nout < 1,
   error('Wrong number of output arguments')
end   
sys = ss(A,B,C,D,1);
%
% end slmoesp
