function [errs,vafs,rcnds] = find_models(y,u,R,list_n,withK)
%FIND_MODELS  Computes the relative error 1-norms and plots pairwisely
%             the estimated output and given output trajectories, for
%             a specified set of system orders, LIST_n.
%
%             FIND_MODELS(Y,U,R,LIST_n)  plots pairwisely the trajectories
%             Y(:,k), and Ye(:,k), for k = 1:l, l = size(Y,2), where Ye is 
%             computed using the given U and estimated system models found
%             by slmoen4 for the set of system orders specified by LIST_n.
%             The argument R is the processed upper triangular factor  R
%             of the concatenated block-Hankel matrices built from the
%             input-output data. R has dimension 2*(m+l)*s-by-2*(m+l)*s,
%             with m = size(U,2). If length(LIST_n) = 2 and LIST_n(2) <= 0,
%             then LIST_n = [ LIST_n(1) : s - 1 ] is used by default.
%
%             [ERRs,VAFs] = FIND_MODELS(Y,U,R,LIST_n)  also returns the
%             relative error 1-norms and the Variance-Accounted-For (VAFs),
%             in percentages, for all outputs of the systems, and all
%             specified orders. ERRs and VAFs have dimensions length(LIST_n)
%             and length(LIST_n)-by-l, respectively.
%
%             [ERRs,VAFs] = FIND_MODELS(Y,U,R,LIST_n,WithK)  also computes
%             and uses the Kalman gain matrix K.
%
%             [ERRs,VAFs,RCNDs] = FIND_MODELS(Y,U,R,LIST_n,WithK)  also
%             returns the reciprocals condition numbers. RCNDs has the
%             dimension length(LIST_n)-by-14, and RCNDs(:,14) contains the
%             forward error bounds for the Riccati equation solutions.
%
%             When there is at least one output argument, plotting the
%             trajectories in the loop list_n can be suppressed by setting
%             the global variable no_loop_plot to 1.
%
%             If requested on output, ERRs and VAFs are always plotted.
%

%        RELEASE 2.0 of SLICOT System Identification Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima 16-09-2000.
%
%        Revisions:
%        V. Sima 30-12-2000, July 2001, March 2009.
%
global pause_wait    % This could be used in pause(n) command.
global no_loop_plot  % Set no_loop_plot = 1 to suppress plotting trajectories.
%
if ~exist('pause_wait', 'var') || isempty(pause_wait),
   pause_wait = -1;  % Standard command pause is used by default.
end
%
if ~exist('no_loop_plot', 'var') || isempty(no_loop_plot),
   no_loop_plot = 0;  % The trajectories are plotted by default.
end
%
nin = nargin;  nout = nargout;
%
if nin < 4,
   disp('Usage: [ERRs,VAFs]       = FIND_MODELS(Y,U,R,LIST_n)')
   disp('       [ERRs,VAFs,RCNDs] = FIND_MODELS(Y,U,R,LIST_n,WithK)')
   return
end
%
l = size(y,2);  m = size(u,2);  s = size(R,1)/(2*( m + l ));
%
if s < 2,
   error('The matrices Y, U, and R have incompatible or too small dimensions')
end
n = list_n(1);  if n < 1,  n = 1;  end
%
if length(list_n) == 2 && list_n(2) <= 0,  list_n = n : s - 1;  end
%
if nin == 4,  withK = 0;  end
%
if l > 4,
   nrplots = fix(sqrt(l));
   plots = [ min( nrplots, 4), min( nrplots, 2) ];
else
   plots = [ min(l,2), 1 ];
end
% 
k = length(list_n);
if nout >= 1,  errs  = zeros(k,1);   end
if nout >= 2,  vafs  = zeros(k,l);   end
if nout == 3,  rcnds = zeros(k,12);  end
%
k = 0;
%
% Loop for the set of orders.
%
for n = list_n,
   k = k + 1;
   if withK == 1,
      [sys,K,rcnd] = slmoen4(s,y,u,n,R);
      [err,ye]     = find_err(y,u,sys,K);
      if nout == 3,  rcnds(k,:) = rcnd';  end
   else
      sys      = slmoen4(s,y,u,n,R);
      [err,ye] = find_err(y,u,sys);
   end
   if nout >= 1,  errs(k) = err;  end
   if nout >= 2,  vafs(k,:) = vaf(y, ye)';  end
   %
   if any( abs(eig(sys.a)) >= 1 ),
      disp(' ')
      disp(['System order n = ', num2str(n),'.  Unstable system !!!'])
      pause(0),
   else
      disp(' ')
      disp(['System order n = ', num2str(n),'.'])
      pause(0),
   end
   %
   if ~no_loop_plot,  plot_ye(y,ye,plots),  end
end
%
if exist('errs', 'var'),
   axis on
   bar(list_n,errs)
   title('Relative estimated output error 1-norms')
   xlabel('System order, n');  ylabel('Error norm')
   if pause_wait < 0,  disp(' ');  disp('Press any key to continue'),  end
   shg,  if pause_wait < 0,  pause,  else  pause(pause_wait),  end 
   close(gcf)
end
%
if exist('vafs', 'var'),
   axis on
   bar(list_n,vafs)
   title('Variance-Accounted-For (all outputs)')
   xlabel('System order, n');  ylabel('VAF')
   if ~exist('errs', 'var') && pause_wait < 0,
      disp(' ');  disp('Press any key to continue'),
   end
   shg,  if pause_wait < 0,  pause,  else  pause(pause_wait),  end  
   close(gcf)
end
%
% end find_models
