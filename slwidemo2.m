echo on,
%   This demonstration example uses a 4-th order model for the flexible
%   robot arm, also identified by slwidemo1.m demonstration file.  
%   Either Cholesky factorization or a Conjugate Gradients algorithm 
%   is used for solving linear systems of equations.
 
echo off
%   RELEASE 2.0 of SLICOT System Identification Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   V. Sima 30-03-2002.
%
%   Revisions: 04-03-2009.
%

close all;

echo on
global pause_wait  % This could be used in pause(n) command.
                   % If pause_wait < 0, standard command pause is used (default).
                   % Any key should then be pressed to continue.

if ~exist('pause_wait', 'var') || isempty(pause_wait),  pause_wait = -1;  end

%       Load the input-output data:

load robot_arm_dat;  u = robot_arm_dat(:,1);   y = robot_arm_dat(:,2);

%       Find a 4-order model based on all data using slmoen4 (alg = 2),
%       using the estimation data set.

est_set = 1 : size(y,1)/2;  val_set = max(est_set)+1 : size(y,1);

alg = 2;  s = 20;  n = 4;  sys = slmoen4(s,y(est_set),u(est_set),n,alg);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

% Estimate the parameters of a Wiener system, using a neural network
% with nn = 12 neurons to model the nonlinear part.
% Default options (except for nprint) and the Cholesky-based 
% Levenberg-Marquardt algorithm are used.
% The error norm is printed in the file IB03AD.prn

nn = 12;  nprint = 1;

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off

ur = u;  yr = y;
ns = max( size( est_set ) );
u  = u(1:ns, :);
y  = y(1:ns, :);

dex_widentc
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off

if nprint == 1 && exist('IB03AD.prn', 'file') == 2,

   % Plot the error norm

   errnrm = textread('IB03AD.prn','%*s%*s%*s%*s%*s%f',-1);
   figure
   set(axes,'FontSize',14)
   plot( errnrm, 'b' );
   title( 'Error norm trajectory in the optimization process' )
end

echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

% Find a possibly better estimate of the nonlinear part,
% setting a tighter tolerance TOL1, and job = 2.

initx = 1;  nprint = 0;  
tol1 = 10^-6;  job = 2;

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off

dex_widentc
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

% Find a possibly better estimate of the whole Wiener system,
% setting a tighter tolerance TOL2, and job = 4.
% The system is initialized by the optimal parameters found so far.

initx = 0;  x = xopt;
tol2 = 10^-6;  job = 4;

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off

dex_widentc
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off

% Reset the main options

initx = [];  job  = [];
tol1  = [];  tol2 = [];

% Now, use default options and the Conjugate Gradients-based 
% Levenberg-Marquardt algorithm.

ialg = 2;

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off

dex_widentc
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off

% Reset ialg

ialg = [];


