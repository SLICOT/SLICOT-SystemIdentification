function theta = s2o(A,B,C,D,x0,apply)
%S2O     Transforms a state-space representation of a linear
%        discrete-time system into the output normal form.
%
%        THETA = S2O(SYS,X0,APPLY)  computes the output normal form,
%        represented by the parameter vector THETA, for a given
%        discrete-time state-space representation of order N, 
%        SYS = (A,B,C,D) (an ss object), with M inputs and L outputs,  
%        with the initial state X0. The vector THETA has the 
%        length N*(L+M+1)+L*M.
%        Instead of the first input parameter SYS, equivalent information
%        may be specified using matrix parameters, for instance,
%        THETA = S2O(A,B,C,D,XO,APPLY);  % or  
%        THETA = S2O(A,B,C,XO,APPLY);    % when D = zeros(L,M);
%
%        THETA = S2O(SYS,APPLY)  assumes that X0 = 0.
%
%        APPLY is an optional integer specifying whether or not the
%        parameter vector THETA should be transformed using a 
%        bijective mapping:
%        APPLY = 1: apply the bijective mapping to the N vectors in
%                   THETA corresponding to the matrices A and C;
%        APPLY = 0: do not apply the bijective mapping.
%        Default:  APPLY = 0.
%        The value of APPLY should be the same as that used in the
%        corresponding call of the paired O2S function.
%
%        See also O2S
%

%        RELEASE 2.0 of SLICOT System Identification Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, Feb. 2002.
%
%        Revisions: Mar. 2009.
%

ni = nargin;
if ni < 1,
   error( 'S2O needs at least 1 input parameters' )
end
if isa( A, 'lti' ),
   % Get the system matrices of the ss object, and the remaining parameters.
   % General call     theta = s2o(A,B,C,D,x0,apply);
   % Special call     theta = s2o(sys,x0,apply); 
   %
   if A.Ts == 0,
      error( 'The system SYS must be a discrete-time system' )
   end
   [ As, Bs, Cs, Ds ] = ssdata( A );
   n = size( As, 1 );
   if ni == 1,
      x0 = zeros( n, 1 );
      theta = ss2onf( As, Bs, Cs, Ds, x0 );
   elseif ni == 2,
      [ nr, nc ] = size( B );  
      if nr*nc  == n,
         % Special call     theta = s2o(sys,x0);  Below, B is x0 ! 
         theta = ss2onf( As, Bs, Cs, Ds, B );
      else
         x0 = zeros( n, 1 );
         apply = B;
         if ~( numel( apply ) == 1 ),
            error( 'APPLY must be an integer variable' )
         end
         if ~( apply == 0 || apply == 1 ),
            apply = 0;  
         end
         theta = ss2onf( As, Bs, Cs, Ds, x0, apply );       
      end
   elseif ni == 3,
      [ nr, nc ] = size( B );  
      if ~( nr*nc  == n ),
         error( 'X0 has wrong size' )
      end
      apply = C;
      if ~( numel( apply ) == 1 ),
         error( 'APPLY must be an integer variable' )
      end
      if ~( apply == 0 || apply == 1 ),
         apply = 0;  
      end
      % Special call     theta = s2o(sys,x0,apply);  Below, B is x0 !
      theta = ss2onf( As, Bs, Cs, Ds, B, apply );
   else
      error( 'Wrong number of input arguments' )
   end
   %
else
   % The system matrices are directly specified.
   % General call     x0 = s2o(A,B,C,D,x0,apply);
   % Special call     x0 = s2o(A,B,C,x0,apply); 
   %                  x0 = s2o(A,B,C,apply); 
   % 
   if ni < 3,
      error( 'S2O needs at least 3 input parameters' )
   end
   [ m2, n2 ] = size( B );  [ m3, n3 ] = size( C );   m = n2;   l = m3;   n = n3;
   if ni >= 4,
      [ m4, n4 ] = size( D );
      if ni >= 5,
         [ m5, n5 ] = size( x0 );
         if ni >= 6,
            if ~( numel( apply ) == 1 ),
                error( 'APPLY must be an integer variable' )
            end
            if ~( apply == 0 || apply == 1 ),
                apply = 0;  
            end
            theta = ss2onf( A, B, C, D, x0, apply );
         else
            % Special calls    theta = s2o(A,B,C,D,x0);
            %                  theta = s2o(A,B,C,x0,apply);
            if m5*n5 == n && m4 == l && n4 == m,             
               theta = ss2onf( A, B, C, D, x0 );
            else             
               % Below, D means x0 !
               theta = ss2onf( A, B, C, zeros( l, m ), D );
            end
         end
      else
         % Special calls    theta = s2o(A,B,C,x0);  
         %                  theta = s2o(A,B,C,apply);
         %                  theta = s2o(A,B,C,D);
         if m4*n4 == n,             
            % Below, D means x0 !
            theta = ss2onf( A, B, C, zeros( l, m ), D );
         elseif m4*n4 == 1,             
            apply = D;
            if ~( apply == 0 || apply == 1 ),
                apply = 0;  
            end
            theta = ss2onf( A, B, C, zeros( l, m ), zeros( n, 1 ), apply );
         elseif m4 == l && n4 == m,             
            x0 = zeros( n, 1 );
            theta = ss2onf( A, B, C, D, x0 );
         else
            error( 'The fourth input argument is wrong' )
         end
      end
   else
      D  = zeros( l, m );
      x0 = zeros( n, 1 );
      theta = ss2onf( A, B, C, D, x0 );
   end
end
%
% end s2o
