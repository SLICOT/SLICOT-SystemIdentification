echo on,
%   This demonstration example uses a 4-th order model for the flexible
%   robot arm, also identified by slidemo1.m demonstration file.  Given
%   this model (the "true system"), output trajectories are generated,
%   by simulation.  Then, based on pairs of input and output trajectories,
%   both deterministic and stochastic systems are identified using
%   SLICOT system identification toolbox, and the results are assessed.
 
echo off
%   RELEASE 2.0 of SLICOT System Identification Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   V. Sima 31-12-2000.
%
%   Revisions: 04-03-2009.
%
echo on
global pause_wait  % This could be used in pause(n) command.
                   % If pause_wait < 0, standard command pause is used (default).
                   % Any key should then be pressed to continue.

if ~exist('pause_wait', 'var') || isempty(pause_wait),  pause_wait = -1;  end

echo off

disp('Deterministic system example:')
disp('-----------------------------') 

echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Load the input-output data:

load robot_arm_dat;  u = robot_arm_dat(:,1);   y = robot_arm_dat(:,2);

%       Find a 4-order model based on all data using slmoen4 (alg = 1).

s = 20;  n = 4;  sys = slmoen4(s,y,u,n);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Print the model:

sys

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Generate an output trajectory using the model and a random input:

time_int = size(u,1);
u = rand(time_int,size(u,2));
y = lsim(sys,u);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Estimate a new model using the input and output trajectories:

syse = slmoen4(s,y,u);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end,  close(gcf)

%       Compute the output error:

x0 = inistate(syse,y,u);
ye = lsim(syse,u,[],x0);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Compute the first nMpar Markov parameters of the two models:

nMpar = 100;

mar_par  = markov(sys, nMpar);
mar_pare = markov(syse,nMpar);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off

disp(' ')
disp( 'Relative output error,')
disp(['                  norm(y - ye)/norm(y)                   = ',...
         sprintf('%0.2e',norm(y - ye)/norm(y))])

disp(' ')
disp( 'Relative error in Markov parameters,')
disp(['                  norm(mar_par - mar_pare)/norm(mar_par) = ',...
         sprintf('%0.2e',norm(mar_par - mar_pare)/norm(mar_par))])

echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off

disp(' ')
disp('Stochastic system example:')
disp('--------------------------') 

echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Perturb the output trajectory:

y = y + 10^-2*rand(time_int,size(y,2));

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Generate another model:

[sysn,K] = slmoen4(s,y,u);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end,  close(gcf)

%       Compute the output errors without, and with a Kalman predictor:

[err,ye]   = find_err(y,u,sysn);
[errK,yeK] = find_err(y,u,sysn,K);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Compute the first nMpar Markov parameters of the new model:

mar_parn = markov(sysn,nMpar);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off

disp(' ')
disp( 'Relative output error,')
disp(['  without predictor:  norm(y - ye)/norm(y)               = ',...
             sprintf('%0.2e',norm(y - ye)/norm(y))])
disp(['  with predictor:     norm(y - yeK)/norm(y)              = ',...
             sprintf('%0.2e',norm(y - yeK)/norm(y))])

disp(' ')
disp( 'Relative error in Markov parameters,')
disp(['                  norm(mar_par - mar_parn)/norm(mar_par) = ',...
         sprintf('%0.2e',norm(mar_par - mar_parn)/norm(mar_par))])

echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Plot the estimated output, without predictor:

plot_ye(y,ye);

%       Plot the estimated output, with predictor:

plot_ye(y,yeK);

%       Compute the autocorrelation of the output error ("residuals"):

res = y - ye;
acc = xcorr(res);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

plot(acc);
title('Autocorrelation of the residuals, without Kalman predictor')
xlabel('Lag');  ylabel('Autocorrelation values')

if pause_wait < 0,  pause,  else  pause(pause_wait),  end,  close(gcf)

%       Also, with Kalman predictor:

resK = y - yeK;
accK = xcorr(resK);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

plot(accK);
title('Autocorrelation of the residuals, with Kalman predictor')
xlabel('Lag');  ylabel('Autocorrelation values')

if pause_wait < 0,  pause,  else  pause(pause_wait),  end,  close(gcf)

echo off
