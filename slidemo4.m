echo on,
%   This demonstration example uses the steel subframe flexible structure
%   data, stored as Example [96-013] of the DAISY collection, available from
%   http://www.esat.kuleuven.ac.be/sista/daisy
%   This is a large data set:  8523 samples; 2 inputs, and 28 outputs.
%   The input signals correspond to two shakers at two locations, and
%   the 28 outputs are accelerations provided by accelerometers around
%   the structure, used for measurements.  The sampling period is
%   1/1024 seconds.  Both inputs are white noise forces.
%
%   The calculations for this example took about 20 minutes for Matlab 
%   N4SID code and over 5 hours for Matlab MOESP (simulation-based), on
%   an IBM-PC computer at 500 MHz with LAPACK-based MATLAB.  All SLICOT 
%   codes (with default settings) solved the problem in about one minute.

echo off
%   RELEASE 2.0 of SLICOT System Identification Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   V. Sima 30-12-2000.
%
%   Revisions: 04-03-2009.
%
echo on
global pause_wait  % This could be used in pause(n) command.
                   % If pause_wait < 0, standard command pause is used (default).
                   % Any key should then be pressed to continue.
                   
global no_loop_plot  % Set no_loop_plot = 1 to suppress plotting trajectories
                     % in the model finding loop.

if ~exist('pause_wait', 'var') || isempty(pause_wait),  pause_wait = -1;  end

%       Load the input-output data:

load flexible_structure_dat;
u = flexible_structure_dat(:,1:2);
y = flexible_structure_dat(:,3:30);
u = detrend(u);  y = detrend(y);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Plot the input-output data:

plot_yu(y,u)

%       Detailed plot for the first 200 samples and first 4 outputs:

plot_yu(y(1:200,1:4),u(1:200,:))

%       Find a model based on all data using slmoen4 (default alg = 1).
%       The model order will be chosen after inspecting the singular values.

s = 21;  n = 20;                    % Could be better to use larger s, n.
[sys,K,rcnd,R] = slmoen4(s,y,u);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end,  close(gcf)

%       Print the model:

sys.a

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

sys.b

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

sys.c

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

sys.d

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Assess the model:

[err(1),ye]  = find_err(y,u,sys);
[err(2),yeK] = find_err(y,u,sys,K);

VAFc = vaf(y, ye);
VAFmean = mean(VAFc);
VAFmin  =  min(VAFc);

VAFcK = vaf(y, yeK);
VAFmeanK = mean(VAFcK);
VAFminK  =  min(VAFcK);

echo off
disp(' ')
disp(['Relative error,       norm(y - ye,1)/norm(y,1)  = ',...
      sprintf('%0.2e',err(1))])
disp(['Idem, with predictor, norm(y - yeK,1)/norm(y,1) = ',...
      sprintf('%0.2e',err(2))])
disp(' ')
if size(y,2) == 1,
   disp('Variance-Accounted-For, in percentages')
   disp(' ')
   disp('            Without predictor   With predictor')
   disp(['mean(VAF) =      ',sprintf('%0.4f',VAFmean),...
              '            ',sprintf('%0.4f',VAFmeanK)])
else
   disp('Means and minimum values of Variance-Accounted-For, in percentages')
   disp('      (for all outputs of a system)')
   disp(' ')
   disp('            Without predictor   With predictor')
   disp(['mean(VAF) =      ',sprintf('%0.4f',VAFmean),...
              '            ',sprintf('%0.4f',VAFmeanK)])
   disp(['min(VAF)  =      ',sprintf('%0.4f',VAFmin),...
              '            ',sprintf('%0.4f',VAFminK)])
end
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Plot the estimated output, without predictor:

plot_ye(y,ye);

%       Plot the estimated output, with predictor:

plot_ye(y,yeK);

%       Detailed plot for the first 200 samples and first 4 outputs,
%       without predictor:

plot_ye(y(1:200,1:4),ye(1:200,1:4))

%       Idem, with predictor:

plot_ye(y(1:200,1:4),yeK(1:200,1:4))

%       Now, check various model orders from 4 to s-1 (step 2), first without the
%       Kalman predictor. The trajectories are not plotted.

list_n = 4 : 2 : 20;  withK = 1;  no_plot_save = no_loop_plot;  no_loop_plot = 1;

[errs,VAFs] = find_models(y,u,R,list_n);

%       Also, check various model orders from 4 to s-1 (step 2), with the
%       Kalman predictor. The trajectories are not plotted.

[errsk,VAFsk] = find_models(y,u,R,list_n,withK);

no_loop_plot = no_plot_save;

%       The goodness of fit increases for increasing n.
%       The system is stable for all orders. 
%
%       Plot the relative error 1-norms and the values for
%       Variance-Accounted-For (VAF) for all orders.

plot(list_n,errs)
title('Relative output error 1-norms without Kalman predictor')
xlabel('System orders');  ylabel('Relative errors')

if pause_wait < 0,  pause,  else  pause(pause_wait),  end,  close(gcf)

plot(list_n,errsk)
title('Relative output error 1-norms with Kalman predictor')
xlabel('System orders');  ylabel('Relative errors')

if pause_wait < 0,  pause,  else  pause(pause_wait),  end,  close(gcf)

plot(list_n,VAFs)
title('Variance-Accounted-For (VAF) values without Kalman predictor')
xlabel('System orders');  ylabel('VAF values')

if pause_wait < 0,  pause,  else  pause(pause_wait),  end,  close(gcf)

plot(list_n,VAFsk)
title('Variance-Accounted-For (VAF) values with Kalman predictor')
xlabel('System orders');  ylabel('VAF values')

if pause_wait < 0,  pause,  else  pause(pause_wait),  end,  close(gcf)

%       Also, print the statistics on Variance-Accounted-For (VAF) values.

echo off
disp(' ')
disp('VAF statistics without Kalman predictor, for all outputs and orders')
disp(' ')
if size(y,2) > 1,
   disp(['mean(VAFs)  = ',sprintf('  %0.2f',mean(VAFs))])
   disp(['min(VAFs)   = ',sprintf('  %0.2f', min(VAFs))])
else
   disp(['mean(VAFs)  = ',sprintf('  %0.2f',mean(VAFs))])
end

disp(' ')
disp('VAF statistics with Kalman predictor, for all outputs and orders')
disp(' ')
if size(y,2) > 1,
   disp(['mean(VAFsk) = ',sprintf('  %0.2f',mean(VAFsk))])
   disp(['min(VAFsk)  = ',sprintf('  %0.2f', min(VAFsk))])
else
   disp(['mean(VAFsk) = ',sprintf('  %0.2f',mean(VAFsk))])
end
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end
echo off
