function plot_yu(y,u)
%PLOT_YU   Plots the output and input trajectories of a system.
%
%          PLOT_YU(Y)  plots the given output trajectories Y.
%
%          PLOT_YU(Y,U)  plots the given output and input
%          trajectories Y and U.

%        RELEASE 2.0 of SLICOT System Identification Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima 30-12-2000.
%
%        Revisions:
%        V. Sima 30-03-2002, 03-03-2009.
%

global pause_wait  % This could be used in pause(n) command.
%
if ~exist('pause_wait', 'var') || isempty(pause_wait),
   pause_wait = -1;  % Standard command pause is used by default.
end
nin = nargin;
%
if nin < 1,
   disp('Usage: PLOT_YU(Y)')
   disp('       PLOT_YU(Y,U)')
   return
end
%
l = size(y,2);
if nin == 2,  m = size(u,2);  else  m = 0;  end
%
axis on
%
if l == m,
   for k = 1 : l,
      subplot(2,1,1), plot(y(:,k)),  xlabel('Samples');  ylabel(['Y_{', num2str(k),'}'])
      subplot(2,1,2), plot(u(:,k)),  xlabel('Samples');  ylabel(['U_{', num2str(k),'}'])
      if ( l > 1 && k < l ) || l == 1,
         if k == 1 && pause_wait < 0,  
            disp(' ');  disp('Press any key to continue'),
         end
         if l > 1,
            shg,  if pause_wait < 0,  pause,  else  pause(pause_wait),  end,  clf,
         end
      end
   end
elseif m == 0,
   for k = 1 : l,
      plot(y(:,k)),  xlabel('Samples');  ylabel(['Y_{', num2str(k),'}'])
      if ( l > 1 && k < l ) || l == 1,
         if k == 1 && pause_wait < 0,  
            disp(' ');  disp('Press any key to continue'),
         end
         if l > 1,
            shg,  if pause_wait < 0,  pause,  else  pause(pause_wait),  end,  clf,
         end
      end
   end
else
   for k = 1 : l,
      for ku = 1 : m,
         subplot(2,1,1), plot(y(:,k)),  xlabel('Samples');  ylabel(['Y_{', num2str(k),'}'])
         subplot(2,1,2), plot(u(:,ku)), xlabel('Samples');  ylabel(['U_{', num2str(ku),'}'])
         if ~( k == l && ku == m ),
            if k == 1 && ku == 1 && pause_wait < 0,
               disp(' ');  disp('Press any key to continue'),
            end
            shg,  if pause_wait < 0,  pause,  else  pause(pause_wait),  end,  clf
         end
      end
   end
end
if pause_wait < 0,  pause,  else  pause(pause_wait),  end,  close(gcf)
%
% end plot_yu
