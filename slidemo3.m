echo on,
%   This demonstration example uses the cutaneous potential recordings
%   of a pregnant woman, stored as Example [96-012] of the DAISY collection,
%   available from http://www.esat.kuleuven.ac.be/sista/daisy
%   There are no inputs, but 8 outputs.  The outputs 1 : 5 are abdominal
%   measurements, and the outputs 6 : 8 are thoracic measurements.
%   The sampling period is 5 seconds.  The number of data samples is 2500.
%   It is highly recommended to use a Kalman predictor with this example.

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

load foetal_ecg_dat;
y = foetal_ecg_dat(:,2:9);
u = zeros(size(y,1),0);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Plot the input-output data:

plot_yu(y)

%       Detailed plot for the first 250 samples:

plot_yu(y(1:250,:))

%       Find a model based on all data using slmoen4 (default alg = 1).
%       The model order will be chosen after inspecting the singular values.

s = 21;
[sys,K,rcnd,R] = slmoen4(s,y);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end,  close(gcf)

%       Print the model:

sys

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

%       Now, check various model orders from 4 to s-1, first without the
%       Kalman predictor. The trajectories are not plotted.

list_n = [4 0];  withK = 1;  no_plot_save = no_loop_plot;  no_loop_plot = 1;

[errs,VAFs] = find_models(y,u,R,list_n);

%       Also, check various model orders from 4 to s-1, with the
%       Kalman predictor. The trajectories are not plotted.

[errsk,VAFsk] = find_models(y,u,R,list_n,withK);

no_loop_plot = no_plot_save;

%       The best fit is obtained for n in the range 12 to 15.
%       The system is unstable for all orders less than 21, except for n = 12
%       and n = 15. 
%
%       Plot the relative error 1-norms and the values for
%       Variance-Accounted-For (VAF) for all orders.

list_n = list_n(1) : s-1;
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

%       The worse identified outputs are 1 and 4, as seen below.

min1 = min(VAFsk(:,1));  ordmin1 = find( min1 == VAFsk(:,1) );
min4 = min(VAFsk(:,4));  ordmin4 = find( min4 == VAFsk(:,4) );
plot(list_n,VAFsk)
title('Variance-Accounted-For (VAF) values with Kalman predictor')
xlabel('System orders');  ylabel('VAF values')
text(ordmin1 + list_n(1) - 1,VAFsk(ordmin1,1),'\leftarrowfor yeK_1')
text(ordmin4 + list_n(1) - 1,VAFsk(ordmin4,4),'\leftarrowfor yeK_4')

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

%       Compute the system model for n = 12, which is stable and provides
%       acceptable good relative errors and VAF values.

n = 12;
[sys,K,rcnd,R] = slmoen4(s,y,u,n);

[err(1),ye]  = find_err(y,u,sys);
[err(2),yeK] = find_err(y,u,sys,K);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Plot the estimated output for n = 12, without predictor:

plot_ye(y,ye);

%       Detailed plot for the first 250 samples:

plot_ye(y(1:250,:),ye(1:250,:))

%       Plot the estimated output for n = 12, with predictor:

plot_ye(y,yeK);

%       Detailed plot for the first 250 samples:

plot_ye(y(1:250,:),yeK(1:250,:))

echo off
