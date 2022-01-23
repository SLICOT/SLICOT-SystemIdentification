% File dex_wident.m
%
% Script using the SLICOT mexfile wident on a given example,
% for estimating the parameters of a Wiener system. 
% The linear part is estimated using SLICOT m-file slmoen4.
% If the option job of wident is specified as 1, 2, or 4, and the 
% variable initx is set to 1, then this m-function computes the
% appropriate part of the initial parameter vector x, but if initx
% is 0, then x should be initialized before calling this script.
% (It is possible to use the results theta and xo of a previous
% initialization with inintx = 1 and the same value of job.)
% The results computed for various values of job will not be always
% identical, since the computations inside and outside wident 
% could slightly differ.
%
% Variables which should be set before calling this script:
% s   - number of block rows; 
% n   - order of the linear part;
% nn  - number of neurons to use;
% ns  - number of samples of the estimation set;
% ur  - the input  trajectory (t-by-m);
% yr  - the output trajectory (t-by-l);
% u   - the part of input  trajectory used for estimation (ns-by-m);
% y   - the part of output trajectory used for estimation (ns-by-l);
% x   - the needed system parameters, if job = 1, 2, or 4, and 
%       initx = 0;
% sys - the Matlab system object for the linear part, if initx = 0.
%
% Variables which could be set before calling this script:
% job    - execution option, for the computation to do (default 3);
% initx  - external initialization option (default 0, i.e., external);
% ITMAX1 - number of of iterations for the initialization of the static
%          nonlinearity (default 500);
% ITMAX2 - number of of iterations for the whole optimization process
%          (default 1000);
% nprint - option specifying the frequency of printing (default 0);
% tol1   - if job = 2 or 3, tolerance for the initialization of the
%          static nonlinearity (default 10^-4);
% tol2   - tolerance for the whole optimization process (default 10^-4);
% seed   - the seed for initializing the random number generator 
%          (default []);
% printw - option for printing warnings (default 1, i.e., print);
% alg    - option for the subspace algorithm for estimating the
%          linear part (default 2, i.e., fast QR);
% range  - the size of moving window for plotting the estimation
%          error trajectories (default 40).

% RELEASE 2.0 of SLICOT System Identification Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% V. Sima, Research Institute for Informatics, Bucharest, Mar. 2002.
%
% Revisions:
% V. Sima, May 2005, Nov. 2005, Mar. 2009.
%
 
close all;

if ~( exist( 'job', 'var' ) ) || isempty( job ),
   job = 3;
end   

if ~( exist( 'initx', 'var' ) ) || isempty( initx ),
   initx = 0;
end   

if ~( exist( 'ITMAX1', 'var' ) ) || isempty( ITMAX1 ),
   ITMAX1 = 500;
end   
if ~( exist( 'ITMAX2', 'var' ) ) || isempty( ITMAX2 ),
   ITMAX2 = 1000;
end   

if ~( exist( 'nprint', 'var' ) ) || isempty( nprint ),
   nprint = 0;
end   

if ~( exist( 'tol1', 'var' ) ) || isempty( tol1 ),
   tol1 = 10^-4;
end   
if ~( exist( 'tol2', 'var' ) ) || isempty( tol2 ),
   tol2 = 10^-4;
end   

if ~( exist( 'seed', 'var' ) ),
   seed = [];
end   

if ~( exist( 'printw', 'var' ) ) || isempty( printw ),
   printw = 1;
end   

if ~( exist( 'alg', 'var' ) ) || isempty( alg ),
   alg = 2;
end   

if exist( 'range', 'var' ) ~= 1,
   range = 40;
elseif isempty( range ),
   range = 40;
end   

[t, m ] = size( ur );
l = size( y, 2 );

if ~exist('pause_wait', 'var') || isempty(pause_wait),  pause_wait = -1;  end

disp( ' ' )
disp( 'Running wident.  Please wait.' )
disp( ' ' )
disp( 'Tolerances for initialization and the whole optimization : ' )
disp( [ 'TOL1 = ', num2str( tol1 ),'  ', 'TOL2 = ', num2str( tol2 ) ] )
disp( ' ' )
disp( 'Maximum number of iterations for initialization and the whole optimization : ' )
disp( [ 'ITMAX1 = ', num2str( ITMAX1 ),'  ', 'ITMAX2 = ', num2str( ITMAX2 ) ] )
disp( ' ')
disp( [ 'Execution option:  job = ',        num2str( job ), '  ', ...
        'Initialization option:  initx = ', num2str( initx ) ] )
disp( ' ')

if ~( job == 3 ) && initx,
    
   % Estimate the linear part and compute the output normal form

   sys = slmoen4( s, y, u, n, alg );
   [A, B, C, D] = ssdata( sys );
   x0 = inistate( sys, y, u );
   ze = ldsim( A, B, C, D, ur, x0 );

   apply = 1;
   theta = ss2onf( A, B, C, D, x0, apply );
   
   % Set known parameters, if any, or compute them, if required
   
   if job == 2 || job == 4,
      disp( 'Setting the linear part parameters' )
      disp( ' ' )
      if job == 2,  
         x = theta;   
      else
         l1 = ( nn*( l + 2 ) + 1 )*l;
         l2 = n*( l + m + 1 ) + l*m;
         x(l1+1:l1+l2) = theta;
      end
   end
   %
   if job == 1 || job == 4,
      disp( 'Initialization of the nonlinear part' )
      disp( ' ' )
      [ xo, perfo, nfo ] = wident( 2, u, y, nn, n, theta, [ITMAX1, 0], ...
                                   nprint, [tol1, tol2], seed, printw );   
      x = xo;
      if printw,
          disp( ' ' )
          disp( ' ' )
      end
   end
else
    
   % Compute the estimated output of the linear part

   if job == 3,  sys = slmoen4( s, y, u, n, alg );  end
   [A, B, C, D] = ssdata( sys );
   x0 = inistate( sys, y, u );
   ze = ldsim( A, B, C, D, ur, x0 );
end

% Optimize the parameters of the whole system

disp( 'Optimization of the whole system' )
disp( ' ' )

tic
   time = cputime;
   if job == 1,
      [ xopt, perf, nf, rcnd ] = wident( job, u, y, nn, s, n, x, [ITMAX1, ITMAX2], ...
                                         nprint, [tol1, tol2], printw );   
   elseif job == 2,
      [ xopt, perf, nf ] = wident( job, u, y, nn, n, x, [ITMAX1, ITMAX2], ...
                                   nprint, [tol1, tol2], seed, printw );   
   elseif job == 3,
      [ xopt, perf, nf, rcnd ] = wident( job, u, y, nn, s, n,  [], [ITMAX1, ITMAX2], ...
                                         nprint, [tol1, tol2], seed, printw );   
   else
      [ xopt, perf, nf ] = wident( job, u, y, nn, n, x, [ITMAX1, ITMAX2], ...
                                   nprint, [tol1, tol2], printw );   
   end
   timings = cputime - time;
   if printw,
       disp( ' ' )
   end
toc

ye = Wiener( n, l, nn, xopt, ur );

errL = zeros(t,1);  errN = zeros(t,1);  errLm = zeros(t,1);  errNm = zeros(t,1);

for k = 1 : t
   errL(k) = norm( yr(k,:) - ze(k,:) );  % error after linear identification
   errN(k) = norm( yr(k,:) - ye(k,:) );  % error after the whole optimization
end;
figure
set(axes,'FontSize',14)
plot( errL, 'r' );
hold
plot( errN, 'g' );
title( 'Errors after linear (red) and Wiener (green) identification' )

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

% Mean value of the error in a moving window of range samples
errLm(1) = sum( errL(1:range) )/range;   
errNm(1) = sum( errN(1:range) )/range;
for k = 1 : t - range
   errLm(k+1) = errLm(k) + ( errL(k + range) - errL(k) )/range;   
   errNm(k+1) = errNm(k) + ( errN(k + range) - errN(k) )/range;   
end;
figure
set(axes,'FontSize',12)
plot( errLm(1:t - range + 1), 'r' );
hold;
plot( errNm(1:t - range + 1), 'g' );
title(  'Mean value of errors for linear and Wiener identification' )
legend( 'Upper curve: linear identification error', ...
        'Lower curve: Wiener identification error' )
xlabel( 'Samples' )

format short e

disp( ' ' )
disp( 'Global performances' )
disp( ' ' )
if job == 2 || job == 3,
   disp( '                          Whole optimization   Initialization of nonlinear part' )
   disp( ' ' )
   disp( ['Sum of squares         :   ', sprintf( '%11.5g', perf(2) ),...
          '                   ',         sprintf( '%11.5g', perf(6) ) ] )
   disp( ['Number of iterations   :   ', sprintf( '%11i',   perf(3) ),...
          '                   ',         sprintf( '%11i',   perf(7) ) ] )
   disp( ['Final Levenberg factor :   ', sprintf( '%11.5g', perf(4) ),...
          '                   ',         sprintf( '%11.5g', perf(8) ) ] )
else
   disp( '                          Whole optimization')
   disp( ' ' )
   disp( [ 'Sum of squares         :   ', sprintf( '%11.5g', perf(2) ) ] )
   disp( [ 'Number of iterations   :   ', sprintf( '%11i',   perf(3) ) ] )
   disp( [ 'Final Levenberg factor :   ', sprintf( '%11.5g', perf(4) ) ] )
end
disp( ' ' )
   
disp( [ 'Total number of function evaluations :   ', sprintf( '%11i', nf(1) ) ] )
disp( [ 'Total number of Jacobian evaluations :   ', sprintf( '%11i', nf(2) ) ] )
disp( ' ')

disp( [ 'Euclidean norm of the error using a linear model :   ', ...
                                            sprintf( '%11.5g', norm( errL ) ) ] )
disp( [ 'Euclidean norm of the error using a Wiener model :   ', ... 
                                            sprintf( '%11.5g', norm( errN ) ) ] )
disp( ' ' )

