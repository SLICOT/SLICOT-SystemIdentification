echo on,
%   This demonstration example uses the heat exchanger example,
%   stored in the DAISY collection,
%   available from http://www.esat.kuleuven.ac.be/sista/daisy
%   There are one input and one output.  The number of data 
%   samples is 4000.

echo off
%   RELEASE 2.0 of SLICOT System Identification Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   V. Sima 30-03-2002.
%
%   Revisions: 04-03-2009.
%

close all

echo on
global pause_wait  % This could be used in pause(n) command.
                   % If pause_wait < 0, standard command pause is used (default).
                   % Any key should then be pressed to continue.
                   
global no_loop_plot  % Set no_loop_plot = 1 to suppress plotting trajectories
                     % in the model finding loop.

if ~exist('pause_wait', 'var') || isempty(pause_wait),  pause_wait = -1;  end

%       Load the input-output data:

load exchanger_dat;
u = exchanger_dat(:,2);  y = exchanger_dat(:,3);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Plot the input-output data:

plot_yu(y,u)

%       Detailed plot for the first 1000 samples:

plot_yu(y(1:1000,:),u(1:1000,:))

%       Find a model based on all data using slmoen4 (alg = 2).
%       The order 6 is used.

alg = 2;  s = 15;  n = 6;  [sys,K,rcnd,R] = slmoen4(s,y,u,n,alg);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Assess the model:

[err(1),ye]  = find_err(y,u,sys);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end,  close(gcf)

%       Plot the estimated output, without predictor:

plot_ye(y,ye);

% Estimate the parameters of a Wiener system, using a neural network
% with nn = 12 neurons to model the nonlinear part.
% Default options (except for nprint) and the MINPACK-like, 
% structure-exploting Levenberg-Marquardt algorithm are used.
% The error norm is printed in the file IB03BD.prn

nn = 12;  nprint = 1;

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off

est_set = 1 : size(y,1)/2;
ur = u;  yr = y;
ns = max( size( est_set ) );
u  = u(1:ns, :);
y  = y(1:ns, :);

dex_wident
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off

if nprint == 1 && exist('IB03BD.prn', 'file') == 2,

   % Plot the error norm

   errnrm = textread('IB03BD.prn','%*s%*s%*s%*s%*s%f',-1);
   figure
   set(axes,'FontSize',14)
   plot( errnrm, 'b' );
   title( 'Error norm trajectory in the optimization process' )
end

% Now, use the Cholesky-based Levenberg-Marquardt algorithm. 
% The error norm is printed in the file IB03AD.prn

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off

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
