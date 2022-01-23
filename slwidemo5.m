echo on,
%   This demonstration example uses a simulated Wiener system,
%   with 3 inputs and 2 outputs.

echo off
%   RELEASE 2.0 of SLICOT System Identification Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   V. Sima 30-03-2002.
%
%   Revisions:
%   V. Sima, Nov. 2005, Mar. 2009.
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

load ex_wident;
ns = 1000;  u = un(1:ns, :);  y = yn(1:ns, :);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Plot the input-output data:

plot_yu(y,u)

%       Detailed plot for the first 250 samples:

plot_yu(y(1:250,:),u(1:250,:))

%       Find a model based on all data using slmoen4 (alg = 2).
%       The model order will be chosen after inspecting the singular values.
%       The order 3 or 4 is suggested to be used.

alg = 2;  s = 10;  [sys,K,rcnd,R] = slmoen4(s,y,u,[],alg);  n = size(sys.a,1);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Assess the model:
  
[err(1),ye] = find_err(y,u,sys);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end,  close(gcf)

%       Plot the estimated output, without predictor:

plot_ye(y,ye);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end,  close(gcf)

%       Also check various orders and plot the estimation errors and VAFs.
%       The VAFs plot indicates that the linear model is not very good,
%       and especially the first output is badly estimated.

list_n = 1 : s-3;
no_loop_plot = 1;  % Set no_loop_plot = 1 to suppress plotting trajectories

[errs,VAFs] = find_models(y,u,R,list_n);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end,  close(gcf)

% Estimate the parameters of a Wiener system, using a neural network
% with nn = 12 neurons to model the nonlinear part.
% Default options (except for nprint) and the MINPACK-like, 
% structure-exploting Levenberg-Marquardt algorithm are used.
% The error norm is printed in the file IB03BD.prn

nn = 12;  nprint = 1;

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off

dex_wident
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off

if nprint == 1 && exist('IB03BD.prn', 'file') == 2,

   % Plot the error norm.
   % The error trajectory includes 3 parts: the first two correpond
   % to individual optimization on the nonlinear part corresponding 
   % to each output, and the third part corresponds to the whole
   % system optimization.

   errnrm = textread('IB03BD.prn','%*s%*s%*s%*s%*s%f',-1);
   figure
   set(axes,'FontSize',14)
   plot( errnrm, 'b' );
   title( 'Error norm trajectory in the optimization process' )
end

echo on

% Find the optimal linear system and see that the output error
% and VAF for this system are still not good enough.

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

l1 = (nn*(l+2)+1)*l;
syso = o2s(n,m,l,xopt(l1+1:end),1);
[erro,yeo] = find_err(y,u,syso);
close all;  plot_ye(y,yeo);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

VAFo = vaf(y,yeo);

echo off

disp( ' ' )
disp( [ 'VAF for the optimal linear system : ', sprintf( '%11.5g', VAFo ) ] )

echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

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
