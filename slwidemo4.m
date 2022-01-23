echo on,
%   This demonstration example uses the ethane-ethylene distillation 
%   column example, stored in the DAISY collection,
%   available from http://www.esat.kuleuven.ac.be/sista/daisy
%   There are 5 inputs and 3 outputs.  The number of data 
%   samples is 4 x 90.

echo off
%   RELEASE 2.0 of SLICOT System Identification Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   V. Sima 30-08-2002.
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

load destill_dat
u1 = destill_dat(:,2:6);    y1 = destill_dat(:,22:24);
u2 = destill_dat(:,7:11);   y2 = destill_dat(:,25:27);
u3 = destill_dat(:,12:16);  y3 = destill_dat(:,28:30);
u4 = destill_dat(:,17:21);  y4 = destill_dat(:,31:33);
u  = [ u1; u2; u3; u4 ];    y  = [ y1; y2; y3; y4 ];   % Concatenation.
ns = size(u,1);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Plot the input-output data:

plot_yu(y,u)

%       Detrend the input-output data:

u = detrend(u);  y = detrend(y);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Find a model based on all data using slmoen4 (alg = 2).
%       The order 3 is used.

alg = 2;  s = 5;  n = 3;  [sys,K,rcnd,R] = slmoen4(s,y,u,n,alg);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Assess the model:
  
[err(1),ye] = find_err(y,u,sys);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end,  close(gcf)

%       Plot the estimated output, without predictor:

plot_ye(y,ye);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end,  close(gcf)

% Estimate the parameters of a Wiener system, using a neural network
% with nn = 12 neurons to model the nonlinear part.
% Default options (except for nprint) and the MINPACK-like, 
% structure-exploting Levenberg-Marquardt algorithm are used.
% The error norm is printed in the file IB03BD.prn

nn = 12;  nprint = 1;
ur = u;   yr = y;

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off

dex_wident
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off

if nprint == 1 && exist('IB03BD.prn', 'file') == 2,
   echo on

   % Plot the error norm.
   % The error trajectory includes 4 parts: the first three correpond
   % to individual optimization on the nonlinear part corresponding 
   % to each output, and the fourth part corresponds to the whole
   % system optimization.

   echo off
   errnrm = textread('IB03BD.prn','%*s%*s%*s%*s%*s%f',-1);
   figure
   set(axes,'FontSize',14)
   plot( errnrm, 'b' );
   title( 'Error norm trajectory in the optimization process' )
end

echo on

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
