function [Y,x] = dsim(A,B,C,D,U,x0)
% DSIM   Computes the output response of a linear discrete-time system.
%
%        [Y,x] = dsim(sys,U,x0)  computes the output vector sequence 
%        y(1), y(2),..., y(t) and the final state x of a discrete-time
%        state-space system, SYS = (A,B,C,D) (an ss object), given the
%        input vector sequence u(1), u(2),..., u(t), and the initial
%        state vector x0. U and Y are matrices with t rows (t is the 
%        number of samples), and as many columns as inputs and outputs,
%        respectively.
%
%        [Y,x] = dsim(sys,U)  uses x0 = 0 as initial state.
%
%        [Y,x] = dsim(A,B,C,D,U,x0)  uses the system matrices instead sys.
%        [Y,x] = dsim(A,B,C,U,x0)    assumes that matrix D is zero.
%
%        See also LDSIM
%

%        RELEASE 2.0 of SLICOT System Identification Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima, March 2002.
%
%        Revisions:
%        V. Sima Mar. 2, 2009.
%

ni = nargin;
if ni < 2,
   error( 'DSIM needs at least 2 input parameters' )
end
if isa( A, 'lti' ),
   % Get the system matrices of the ss object, and the remaining parameters.
   % General call     [Y,x] = dsim(A,B,C,D,U,x0);
   % Special call     [Y,x] = dsim(sys,U,x0); 
   %
   if A.Ts == 0,
      error( 'The system SYS must be a discrete-time system' )
   end
   [ As, Bs, Cs, Ds ] = ssdata( A );
   [ l, m ] = size( A );  n = size( As, 1 );
   if ~( size( B, 2 ) == m ),
      error( 'The matrix U must have as many columns as inputs' )
   end
   if ni == 2,  
      % Special call     [Y,x] = dsim(sys,U);  Below, B is U ! 
      [ Y, x ] = ldsim( As, Bs, Cs, Ds, B );
   elseif ni == 3,
      if numel( C ) < n,
         error( 'The initial state vector is too short' )
      end
      % Special call     [Y,x] = dsim(sys,U,x0);  Below, B is U and C is x0 ! 
      [ Y, x ] = ldsim( As, Bs, Cs, Ds, B, C );
   else
      error( 'Wrong number of input arguments' )
   end
   %
else
   % The system matrices are directly specified.
   % General call     [Y,x] = dsim(A,B,C,D,U,x0);
   % Special calls    [Y,x] = dsim(A,B,C,U,x0); 
   %                  [Y,x] = dsim(A,B,C,U); 
   % 
   if ni < 4,
      error( 'DSIM needs at least 4 input parameters' )
   end
   [ m2, n2 ] = size( B );  [ m3, n3 ] = size( C );   m = n2;   l = m3;   n = n3;
   [ m4, n4 ] = size( D );
   if ni >= 5,
      [ m5, n5 ] = size( U );
      if ni >= 6,
         if numel( x0 ) < n,
            error( 'The initial state vector is too short' )
         end
         [ Y, x ] = ldsim( A, B, C, D, U, x0 );
      else
         % Special calls    [Y,x] = dsim(A,B,C,D,U);
         %                  [Y,x] = dsim(A,B,C,U,x0);
         if m4 == l && n4 == m && n5 == m,             
            [ Y, x ] = ldsim( A, B, C, D, U );
         else             
            % Below, D means U, and U means x0 !
            [ Y, x ] = ldsim( A, B, C, zeros( l, m ), D, U );
         end
      end
   else
      % Special call     [Y,x] = dsim(A,B,C,U);  
      if n4 == m,
         % Below, D means U !
         [ Y, x ] = ldsim( A, B, C, zeros( l, m ), D );
      else
         error( 'The fourth input argument is wrong' )
      end
   end
end
%
% end dsim
