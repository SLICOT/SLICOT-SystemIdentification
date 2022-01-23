%   The SLICOT Wiener System Identification Toolbox contains tools for
%   finding mathematical models of nonlinear Wiener multivariable
%   discrete-time dynamical systems, using input-output data.
%   The linear part is identified using fast subspace-based techniques 
%   and the state-space model is converted to the output normal form.
%   The nonlinear part is modelled as a single layer neural network.
%   All parameters are optimized using Levenberg-Marquardt algorithms.
%
%   SLICOT Wiener Identification Toolbox demonstrations:
%
%   1) A simple example: a flexible robot arm.
%   2) The flexible robot arm: Cholesky or conjugate gradients solver.
%   3) Heat exchanger.
%   4) Ethane-ethylene distillation column.
%   5) System with 3 inputs and 2 outputs.
%      (It could take about 5 minutes on a 500 MHz machine.)
%
%   0) Quit.
%
%   Note: SLICOT Wiener System Identification Toolbox uses the system object 
%         from Matlab Control System Toolbox.  

%   RELEASE 2.0 of SLICOT System Identification Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   V. Sima 30-03-2002.
%
%   Revisions:
%   V. Sima 03-07-2020.
%

global pause_wait  % This could be used in pause(n) command.
                   % If pause_wait < 0, standard command pause is used (default).
                   % Any key should then be pressed to continue.
                   % It may need to use the command "pause on" before
                   % calling fstdemo.

k = 0;
while (1)
   disp(' ')
   help slwidemo
   k = input('Select a demonstration example number: ');
   disp(' ')
   if isempty(k), k = 20;  end
   if k == 0,  break,      end
   if k == 1,  slwidemo1,  end
   if k == 2,  slwidemo2,  end
   if k == 3,  slwidemo3,  end
   if k == 4,  slwidemo4,  end
   if k == 5,  slwidemo5,  end
   close all
end
