function plot_ye(y,ye,plots)
%PLOT_YE   Plots pairwisely the estimated output and given output
%          trajectories.
%
%          PLOT_YE(Y,Ye)  plots pairwisely the trajectories
%          Y(:,k), and Ye(:,k), for k = 1:l, l = size(Y,2).
%
%          PLOT_YE(Y,Ye,PLOTS)  the figure window is broken into
%          maximum PLOTS(1)-by-PLOTS(2) subplots.
%          Default: PLOTS = [min(2,l),1].
%

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
if nin < 2
   disp('Usage: PLOT_YE(Y,Ye)')
   disp('       PLOT_YU(Y,Ye,PLOTS)')
   return
end
%
l = size(y,2);
if nin < 3 || isempty(plots),  plots = [min(2,l),1];  end
%
axis on
nplots = min( plots(1)*plots(2), l );
nplcrt = nplots;
npl    = mod( l, nplots );
nloops = fix( l/nplots );
if ~( npl == 0 ),
   nloops = nloops + 1;
else
   npl = nplots;
end
%
for k = 1 : nloops,
   if k == nloops,  nplcrt = npl;  end
   for np = 1 : nplcrt,
      j = ( k - 1 )*nplots + np;
      err2 = norm(y(:,j) - ye(:,j))/max( eps, norm(y(:,j)));
      subplot(plots(1),plots(2),np), plot([y(:,j) ye(:,j)])
      if plots(1) <= 2,
         if plots(2) <= 2,
            title(['norm(Y_{',num2str(j),'}-Ye_{',num2str(j),...
                '})/norm(Y_{',num2str(j),'}) = ' ,num2str( err2 )])
         elseif plots(2) <= 4,
            title(['Relerr_{',num2str(j),'} = '  ,num2str( err2 )])
         end
      end
      if ismember( np, ( nplcrt - plots(2) + 1 : nplcrt ) ),  xlabel('Samples');  end
      ylabel(['Y_{', num2str(j), '},', 'Ye_{', num2str(j), '}'])
   end
   if k < nloops || nloops == 1,
      if k == 1 && pause_wait < 0,
         disp(' ');  disp('Press any key to continue'),
      end
      if nloops > 1,
         shg,  if pause_wait < 0,  pause,  else  pause(pause_wait),  end,  clf,
      end
   end
end
if pause_wait < 0,  pause,  else  pause(pause_wait),  end,  close(gcf)
%
% end plot_ye
