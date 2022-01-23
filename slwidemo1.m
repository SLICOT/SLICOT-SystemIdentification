echo on,
%   This demonstration example uses the data for a flexible robot arm,
%   stored as Example [96-009] of the DAISY collection, available from
%   http://www.esat.kuleuven.ac.be/sista/daisy
%   The arm is installed on an electrical motor.  The output y is
%   the measured reaction torque of the structure on the ground.
%   The input u is the acceleration of the flexible arm.  The applied
%   input is a periodic sine sweep.  The number of data samples is 1024.
 
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

load robot_arm_dat
u = robot_arm_dat(:,1);
y = robot_arm_dat(:,2);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Plot the input-output data:

plot_yu(y,u)

%       Detailed plot for the first 200 samples:

plot_yu(y(1:200),u(1:200))

%       Define an estimation data set and a validation data set,
%       as the first half and the last half of the time interval:

est_set = 1 : size(y,1)/2;  val_set = max(est_set)+1 : size(y,1);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Find a model based on the estimation data set using slmoen4 (alg = 2).
%       A 4-order model is good enough.

alg = 2;  s = 20;   trial = 1;
[sys,K,rcnd,R] = slmoen4(s,y(est_set),u(est_set),[],alg);

%       Print the model:

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

sys

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Assess the model on the validation data set:

[err(trial,1),ye] = find_err(y(val_set),u(val_set),sys);

VAFc = vaf(y(val_set), ye);
VAFmean(trial) = mean(VAFc); 
VAFmin( trial) =  min(VAFc); 

echo off
disp(' ')
disp(['Relative error,       norm(y - ye,1)/norm(y,1)  = ',...
      sprintf('%0.2e',err(trial,1))])
disp(' ')
if size(y,2) == 1,
   disp('Variance-Accounted-For, in percentages')
   disp(' ')
   disp(['VAF = ',num2str(VAFmean(trial))])
else
   disp('Means and minimum values of Variance-Accounted-For, in percentages')
   disp('      (for all outputs of a system)')
   disp(['mean(VAF) = ',num2str(VAFmean(trial))])
   disp(['min(VAF)  = ',num2str(VAFmin(trial))])
end
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Plot the estimated output:

plot_ye(y(val_set),ye);

%       Now, check various model orders from 4 to 7. 
%       The 6-th order model is slightly better.

list_n = 4 : 7;

[errs,VAFs] = find_models(y,u,R,list_n);
echo off

%       Select the order of the linear part. 

disp(' ')
n = input('Select the order of the linear part (n = 4 is suggested), n = ');
echo on

% Estimate the parameters of a Wiener system, using a neural network
% with nn = 12 neurons to model the nonlinear part.
% Default options (except for nprint) and the MINPACK-like, 
% structure-exploting Levenberg-Marquardt algorithm are used.
% The error norm is printed in the file IB03BD.prn

nn = 12;  nprint = 1;

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off

ur = u;  yr = y;
ns = max( size( est_set ) );
u  = u(1:ns, :);
y  = y(1:ns, :);

dex_wident
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off

if nprint == 1 && exist('IB03BD.prn', 'file') == 2,
   echo on

   % Plot the error norm

   echo off
   errnrm = textread('IB03BD.prn','%*s%*s%*s%*s%*s%f',-1);
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

dex_wident
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

% Find a possibly better estimate of the whole Wiener system,
% setting a tighter tolerance TOL2, and job = 4.
% The system is initialized by the optimal parameters found so far.

initx = 0;  x = xopt;
tol2 = 10^-6;  job = 4;

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off

dex_wident
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off

% Reset the main options

initx = [];  job  = [];
tol1  = [];  tol2 = [];
