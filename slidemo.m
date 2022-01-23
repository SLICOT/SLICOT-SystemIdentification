%   The SLICOT System Identification Toolbox contains subspace-based
%   tools for finding mathematical models of linear multivariable 
%   dynamical systems, using input-output data.
%
%   SLICOT Identification Toolbox demonstrations:
%
%   1) A simple example: a flexible robot arm.
%   2) Given "true" system: the flexible robot arm model in 1).
%   3) Complex system without inputs: cutaneous potential recordings
%      of a pregnant woman.
%   4) More complex system with 28 outputs: steel subframe flexible structure.
%      (It could take about 5 minutes on a 500 MHz machine.)
%   5) Compare different SLICOT identification methods and algorithms.
%
%   0) Quit.
%
%   Note: SLICOT System Identification Toolbox uses the system object 
%         from Matlab Control System Toolbox.  The demonstration for 2)
%         also uses the function xcorr from Signal Processing Toolbox.

%   RELEASE 2.0 of SLICOT System Identification Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   V. Sima 30-12-2000.
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
   help slidemo
   k = input('Select a demonstration example number: ');
   disp(' ')
   if isempty(k), k = 20; end
   if k == 0,  break,     end
   if k == 1,  slidemo1,  end
   if k == 2,  slidemo2,  end
   if k == 3,  slidemo3,  end
   if k == 4,  slidemo4,  end
   if k == 5,  slidemo5,  end
end
