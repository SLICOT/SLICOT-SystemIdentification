echo on,
%   This demonstration example uses the data for a 120 MW power plant,
%   stored as Example [96-003] of the DAISY collection, available from
%   http://www.esat.kuleuven.ac.be/sista/daisy
%   The process has 3 outputs and 5 inputs.  The outputs are: steam 
%   pressure, main steam temperature, and reheat steam temperature.
%   The inputs are: gas flow, turbine valves opening, super heater spray
%   flow, gas dampers, and air flow.  The number of data samples is 200
%   and the sampling time is 1228.8 seconds.
%   SLICOT function slmoen4 is used first, and then other SLICOT 
%   functions and algorithms are called.
 
echo off
%   RELEASE 2.0 of SLICOT System Identification Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   V. Sima 30-12-2000.
%
%   Revisions: 04-03-2009, 08-01-2018.
%
echo on
global pause_wait  % This could be used in pause(n) command.
                   % If pause_wait < 0, standard command pause is used (default).
                   % Any key should then be pressed to continue.

if ~exist('pause_wait', 'var') || isempty(pause_wait),  pause_wait = -1;  end

%       Load the input-output data:

load powerplant_dat
u = powerplant_dat(:,2:6);
y = powerplant_dat(:,7:9);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Plot the input-output data:

plot_yu(y,u)

%       Find a model using slmoen4 (default alg = 1).
%       A 6-order model is good enough.

s = 10;
[sys,K,rcnd,R] = slmoen4(s,y,u);

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
   disp(['VAF = ',num2str(VAFmean)])
else
   disp('Means and minimum values of Variance-Accounted-For, in percentages')
   disp('      (for all outputs of a system)')
   disp(' ')
   disp(['mean(VAF) = ',num2str(VAFmean)])
   disp(['min(VAF)  = ',num2str(VAFmin)])
end
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Plot the estimated output, without predictor:

plot_ye(y,ye);

%       Plot the estimated output, with predictor:

plot_ye(y,yeK);

%       Now, check various model orders from 4 to 9, first without, and then
%       with the Kalman predictor.
%       The system is unstable for n = 4 and n = 5.

list_n = 4 : 9;  withK = 1;

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

%       Now, check other SLICOT functions and algorithms for n = 6.
%       Modulo a similarity transformation, all computed systems have
%       practically the same matrices A, C and K;  for slmoen4 and
%       sln4sid, this is also true for the matrices B and D.
%       First, compute a reference system, using slmoen4, as before.

n = 6;  max_err = 0;
[sys,K,rcnd,R] = slmoen4(s,y,u,n);
 
nrm_A = max([ eps, norm(sys.a,1)]);  
nrm_B = max([ eps, norm(sys.b,1)]);
nrm_C = max([ eps, norm(sys.c,1)]);  
nrm_D = max([ eps, norm(sys.d,1)]);
nrm_K = max([ eps, norm(K,1)]);

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Loop for methods and algorithms.
%       The relative output errors are also stored.

err  = zeros(3,4);
errK = zeros(3,4);

for alg = 1 : 3,
   if alg == 2,  echo off,  end
   [sys1,K1] = slmoen4(s,y,u,n,alg);
   err( alg,1) = find_err(y,u,sys1);
   errK(alg,1) = find_err(y,u,sys1,K1);
   nrm_err = max([ norm(abs(sys1.a) - abs(sys.a),1)/nrm_A, ... 
		   norm(abs(sys1.b) - abs(sys.b),1)/nrm_B, ...
		   norm(abs(sys1.c) - abs(sys.c),1)/nrm_C, ...
		   norm(abs(sys1.d) - abs(sys.d),1)/nrm_D, ...
		   norm(abs(K1)     - abs(K),1)    /nrm_K]);
   max_err = max( max_err, nrm_err);
   
   [sys2,K2] = sln4sid(s,y,u,n,alg);
   err( alg,2) = find_err(y,u,sys2);
   errK(alg,2) = find_err(y,u,sys2,K2);
   nrm_err = max([ norm(abs(sys2.a) - abs(sys.a),1)/nrm_A, ... 
		   norm(abs(sys2.b) - abs(sys.b),1)/nrm_B, ...
		   norm(abs(sys2.c) - abs(sys.c),1)/nrm_C, ...
		   norm(abs(sys2.d) - abs(sys.d),1)/nrm_D, ...
		   norm(abs(K2)     - abs(K),1)    /nrm_K]);
   max_err = max( max_err, nrm_err);
   
   [sys3,K3] = slmoesp(s,y,u,n,alg);
   err( alg,3) = find_err(y,u,sys3);
   errK(alg,3) = find_err(y,u,sys3,K3);
   nrm_err = max([ norm(abs(sys3.a) - abs(sys.a),1)/nrm_A, ... 
		   norm(abs(sys3.c) - abs(sys.c),1)/nrm_C, ...
		   norm(abs(K3)     - abs(K),1)    /nrm_K]);
   max_err = max( max_err, nrm_err);
   
   [sys4,K4] = slmoesm(s,y,u,n,alg);
   err( alg,4) = find_err(y,u,sys4);
   errK(alg,4) = find_err(y,u,sys4,K4);
   nrm_err = max([ norm(abs(sys4.a) - abs(sys.a),1)/nrm_A, ... 
		   norm(abs(sys4.c) - abs(sys.c),1)/nrm_C, ...
		   norm(abs(K4)     - abs(K),1)    /nrm_K]);
   max_err = max( max_err, nrm_err);
end   

echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Print the maximum relative error in system matrices and Kalman
%       gain, compared to the reference results (got using slmoen4), 
%       for all method-oriented SLICOT functions and algorithms.

echo off
disp(' ')
disp(['Maximum relative error in system matrices and Kalman gain : ',...
       num2str(max_err)])

echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Plot the relative errors for all functions and algorithms.
%       First, without Kalman gain.

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

set(axes,'FontSize',16)
bar(err)
legend('slmoen4','sln4sid','slmoesp','slmoesm','location','Best')
xlabel('Algorithm # (1 Cholesky, 2 Fast QR, 3 QR)')
ylabel('Relative errors (without Kalman gain)')

shg,  if pause_wait < 0,  pause,  else  pause(pause_wait),  end
close(gcf)

%       Then, with Kalman gain.

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

set(axes,'FontSize',16)
bar(errK)
legend('slmoen4','sln4sid','slmoesp','slmoesm','location','Best')
xlabel('Algorithm # (1 Cholesky, 2 Fast QR, 3 QR)')
ylabel('Relative errors (with Kalman gain)')

shg,  if pause_wait < 0,  pause,  else  pause(pause_wait),  end  
close(gcf)

%       Finally, a command for obtaining the state, output, and
%       state-output (cross-)covariance matrices Q, Ry, and S, is shown.
%       A combined method is used: A and C via MOESP, B and D via N4SID.

[N,l] = size(y);  meth = 3;

[sys1,K1,Q,Ry,S,rcnd] = findABCD(s,n,l,R,meth,N);

echo off

Q

echo on
if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off

Ry

echo on
if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off

S

echo on
if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
