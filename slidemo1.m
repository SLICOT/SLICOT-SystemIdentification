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
%   V. Sima 30-12-2000.
%
%   Revisions: 04-03-2009.
%
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

%       Find a model based on all data using slmoen4 (default alg = 1).
%       A 4-order model is good enough.

s = 20;  trial = 1;
[sys,K,rcnd,R] = slmoen4(s,y,u);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end,  close(gcf)

%       Print the model:

sys

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Assess the model:

[err(trial,1),ye]  = find_err(y,u,sys);
[err(trial,2),yeK] = find_err(y,u,sys,K);

VAFc = vaf(y, ye);
VAFmean(trial) = mean(VAFc);
VAFmin( trial) =  min(VAFc);

echo off
disp(' ')
disp(['Relative error,       norm(y - ye,1)/norm(y,1)  = ',...
      sprintf('%0.2e',err(trial,1))])
disp(['Idem, with predictor, norm(y - yeK,1)/norm(y,1) = ',...
      sprintf('%0.2e',err(trial,2))])
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

%       Plot the estimated output, without predictor:

plot_ye(y,ye);

%       Plot the estimated output, with predictor:

plot_ye(y,yeK);

%       Now, define an estimation data set and a validation data set,
%       as the first half and the last half of the time interval:

est_set = 1 : size(y,1)/2;  val_set = max(est_set)+1 : size(y,1);

%       Find a 4-order model based on estimation data set using slmoen4:

trial = 2;  n = 4;
[sys,K,rcnd] = slmoen4(s,y(est_set),u(est_set),n);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Assess the model on the validation data set:

[err(trial,1),yeh]  = find_err(y(val_set),u(val_set),sys);
[err(trial,2),yeKh] = find_err(y(val_set),u(val_set),sys,K);

VAFch = vaf(y(val_set), yeh);
VAFmean(trial) = mean(VAFch); 
VAFmin( trial) =  min(VAFch); 

echo off
disp(' ')
disp(['Relative error,       norm(y - ye,1)/norm(y,1)  = ',...
      sprintf('%0.2e',err(trial,1))])
disp(['Idem, with predictor, norm(y - yeK,1)/norm(y,1) = ',...
      sprintf('%0.2e',err(trial,2))])
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

%       Plot the estimated output, without predictor:

plot_ye(y(val_set),yeh);

%       Plot the estimated output, with predictor:

plot_ye(y(val_set),yeKh);

%       The accuracy of the model does not deteriorate when using
%       the estimation data set.  The rows of the matrix err contain
%       the relative estimated output error 1-norms when using all data
%       and estimation data set, respectively.  The first column
%       contains the error norm values without a Kalman predictor
%       and the second column, the corresponding values with predictor.

echo off
err
echo on

%       Now, check various model orders from 4 to 7, first without, and then
%       with the Kalman predictor. The 6-th order model is slightly better.
%       Larger orders have also been checked.  For n = 12 and with K,
%       a perfect fit is got.  The system is unstable for n = 8 : 10, and
%       n > 15.

list_n = 4 : 7;  withK = 1;

[errs,VAFs] = find_models(y,u,R,list_n);

[errsk,VAFsk] = find_models(y,u,R,list_n,withK);

%       Also, print relative error 1-norms and statistics on
%       Variance-Accounted-For (VAF) values.

echo off
disp(' ')
disp(['Relative output error 1-norms for system orders n = [ ',...
       num2str(list_n(1)),' : ',num2str(list_n(length(list_n))),' ]'])
disp(' ')
disp('Without Kalman predictor')
errs

echo on
if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' ')
disp('With Kalman predictor')
errsk

echo on
if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp(' ')
disp('VAF statistics without Kalman predictor')
disp(' ')
if size(y,2) > 1,
   disp(['Mean VAF value for all outputs       = ', num2str(mean(VAFs))])
   disp(['Mimimum VAF value for all outputs    = ', num2str(min(VAFs))])
else
   disp(['mean(VAF values) = ', num2str(mean(VAFs))])
end

disp(' ')
disp('VAF statistics with Kalman predictor')
disp(' ')
if size(y,2) > 1,
   disp(['Mean VAF value for all outputs       = ', num2str(mean(VAFsk))])
   disp(['Mimimum VAF value for all outputs    = ', num2str(min(VAFsk))])
else
   disp(['mean(VAF values) = ', num2str(mean(VAFsk))])
end

echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
