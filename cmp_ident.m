% System identification performance demo using data from the Daisy collection:
% Comparison between SLICOT function slmoen4 and MATLAB function n4sid.
%
% The SLICOT variables alg, printw, conct, and tol can be specified.
% Default: alg = 1, printw = 0, conct = 2, and tol = [0,-1].
% Performance variables set: timing, err, VAFmean, VAFmin.
%
% For the MATLAB n4sid function, with option 'N4Weight' = 'MOESP', the 
% variables given_n, covr, given_N4H, no, and N4H can be specified.
% Default: given_n = 1, covr = 'none', given_N4H = 0,
%          no = zeros(4,1), N4H = zeros(4,3).
% Note that for the value covr = 'none', the parameter covariances are 
% not computed.  To compute them (could be very slow), set covr = [].
% Performance variables set: timing, err, VAFmean, VAFmin, no, N4H.
%
% Set given_n   = 0 for finding the best order in 1 : 10;
%     given_n   = 1 for using the given order n;
%     given_n   = 2 for using the order set in the corresponding no.
% Set given_N4H < 0 for computing minimum horizons;
%     given_N4H = 0 for computing default horizons ('N4Horizon' = 'Auto', or []);
%     given_N4H = 1 for using [s s s];
%     given_N4H = 2 for using given horizons (stored in N4H).
%

%        RELEASE 2.0 of SLICOT System Identification Toolbox.
%        Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%        V. Sima Nov. 2004.
%
%        Revisions:
%        V. Sima Jan. 2007, Mar. 2009, 08-01-2018.
%

disp( ' ' );
disp( 'System Identification performance demo using some data from the Daisy collection' );
disp( ' ' );
disp( 'Applications:' );
disp( ' 1) Industrial evaporator (6305 samples; 3 inputs, 3 outputs)' );
disp( ' 2) Continuous stirred tank reactor (7500 samples; 1 inputs, 2 outputs)' );
disp( ' 3) Steam generator (9600 samples; 4 inputs, 4 outputs)' );
disp( ' 4) CD-player arm (2048 samples; 2 inputs, 2 outputs)' );
disp( ' ' );
%
nofexampl = 4;                           % Number of examples.
%
dim = zeros(nofexampl,6);

err     = repmat(NaN,nofexampl,2);
VAFmean = repmat(NaN,nofexampl,2);
VAFmin  = repmat(NaN,nofexampl,2);
timing  = repmat(NaN,nofexampl,2);
%
examples = [4 10 11 14];                 % Example list using numbering in (Sima et. al, 2004).
%
% SLICOT options.
%
if ~exist( 'alg', 'var' ) || isempty( alg )
   alg = 1;                              % Cholesky factorization algorithm.
end
if ~exist( 'printw', 'var' ) || isempty( printw )
   printw = 0;                           % Do not print SLICOT warnings.
end
if ~exist( 'conct', 'var' ) || isempty( conct )
   conct = 2;                            % No connection between the data blocks.
end
if ~exist( 'tol', 'var' ) || isempty( tol )
   tol = [0; -1];                        % Default SLICOT tolerances.
end
%
% MATLAB options.
%
if ~exist( 'given_n', 'var' ) || isempty( given_n )
   given_n = 1;                          % Order n is given (same as for SLICOT).
end
if ~exist( 'given_N4H', 'var' ) || isempty( given_N4H )
   given_N4H = 0;                        % Horizons automatically chosen.
end
if ~exist( 'covr', 'var' ) || isempty( covr )
   covr = 'none';                        % Parameter covariances not computed.
end
if ~exist( 'no', 'var' ) || isempty( no ) 
   no = zeros( nofexampl, 2 );           % Set to order used.
end
if ~exist( 'N4H', 'var' ) || isempty( N4H ) 
   N4H = zeros( nofexampl, 3 );          % Set to horizons used.
end
%
disp( 'SLICOT calculations' );
disp( '-------------------' );
disp( ' ' )
disp( ['Algorithm parameter, alg = ', sprintf( '%3d', alg )] )
%
appl = 0;
for example = examples
   sample = example;
   appl = appl + 1;
   disp(' ')
   disp(['Application # ', num2str(appl)])
   %
   switch sample
      case 4,
         load evaporator_dat         % D004  n = 4
         u = evaporator_dat(:,1:3);
         y = evaporator_dat(:,4:6);
         s = 10;  n = 4;
         clear evaporator_dat
      case 10,
         load cstr_dat               % D010  n = 5
         u = cstr_dat(:,2);
         y = cstr_dat(:,3:4);
         u = detrend(u);  y = detrend(y);
         s = 15;  n = 5;
         clear cstr_dat
      case 11,
         load steamgen_dat           % D011  n = 9
         u = steamgen_dat(:,2:5);
         y = steamgen_dat(:,6:9);
         u = detrend(u);  y = detrend(y);
         s = 15;  n = 9;
         clear steamgen_dat
      case 14,
         load CD_player_arm_dat      % D103  n = 8
         u = CD_player_arm_dat(:,1:2);
         y = CD_player_arm_dat(:,3:4);
         u = detrend(u);  y = detrend(y);
         s = 15;  n = 8;
         clear CD_player_arm_dat
   end
   %
   [t,l] = size( y );  m = size( u, 2 );
   dim(appl,:) = [appl t m l s n ];
   if ~( sample == 100 ),
      %
      lc = 1;
      time = cputime;
      [ sys, K ] = slmoen4( s, y, u, n, alg, tol, printw );
      timing(appl,lc) = cputime - time;
      [errc,ye] = find_err(y,u,sys,K);
      err(appl,lc) = errc; 
      VAFc = vaf(y, ye);
      VAFmean(appl,lc) = mean(VAFc); 
      VAFmin(appl,lc)  = min(VAFc); 
      %
   end
   clear K rcnd sys u u1 u2 u3 u4 y y1 y2 y3 y4 ye VAFc
end
%
disp( ' ' );
disp( ' ' );
disp( 'MATLAB calculations' );
disp( '-------------------' );
disp( ' ' )
disp( 'Options for n4sid:' )
if given_n == 0,  
   disp( '  order = ''best''' )
elseif given_n == 1,  
   disp( '  Given order' )
else
   disp( '  Order set in the corresponding no' )
end
disp( '  ''N4Weight'' = ''MOESP''' )
if given_N4H < 0,  
   disp( '  Minimum horizons: N4Horizon = [0,0,0]' )
elseif given_N4H == 0,  
   disp( '  N4Horizon = ''Auto''' )
elseif given_N4H == 1,  
   disp( '  N4Horizon = [s s s]' )
else
   disp( '  N4Horizon stored in N4H' )
end
disp( ['  CovarianceMatrix = ''', covr, ''''] )
%
appl = 0;
for example = examples
   sample = example;
   appl = appl + 1;
   disp(' ')
   disp( ['Application # ', num2str( appl )] )
   %
   switch sample
      case 4,
         load evaporator_dat         % D004  n = 4
         u = evaporator_dat(:,1:3);
         y = evaporator_dat(:,4:6);
         s = 10;  n = 4;
         clear evaporator_dat
      case 10,
         load cstr_dat               % D010  n = 5
         u = cstr_dat(:,2);
         y = cstr_dat(:,3:4);
         u = detrend(u);  y = detrend(y);
         s = 15;  n = 5;
         clear cstr_dat
      case 11,
         load steamgen_dat           % D011  n = 9
         u = steamgen_dat(:,2:5);
         y = steamgen_dat(:,6:9);
         u = detrend(u);  y = detrend(y);
         s = 15;  n = 9;
         clear steamgen_dat
      case 14,
         load CD_player_arm_dat      % D103  n = 8;
         u = CD_player_arm_dat(:,1:2);
         y = CD_player_arm_dat(:,3:4);
         u = detrend(u);  y = detrend(y);
         s = 15;  n = 8;
         clear CD_player_arm_dat
   end
   %
   [t,l] = size( y );  m = size( u, 2 );
   lc = 2;
   %
   % MATLAB options.
   %
   if given_n == 0,  
      ord = 'best';  
   elseif given_n == 1,  
      ord = n;  
   else
      ord = no(appl,lc);  
   end
   if given_N4H < 0,  
      N4Hc = 0;  
   elseif given_N4H == 0,  
      N4Hc = 'Auto';  
   elseif given_N4H == 1,  
      N4Hc = [s s s];  
   else
      N4Hc = N4H(appl,1:3);  
   end
   %
   time = cputime;
   data = iddata( y, u );
   model = n4sid( data, ord, 'n4w', 'moesp', 'N4H', N4Hc, 'cov', covr );
   if size( u, 2 ) == 0
       [A,B,C,D] = ssdata( model );
       sys = ss( A, B, C, D, 1 );
   else
       sys = ss( model( 'm' ) );
   end
   timing(appl,lc) = cputime - time;
   no(appl,lc)     = size( sys.a, 1 );
   N4H(appl,1:3)   = model.EstimationInfo.N4Horizon;
   K = model.k;
   [errc,ye] = find_err( y, u, sys, K );
   err(appl,lc) = errc; 
   VAFc = vaf(y, ye);
   VAFmean(appl,lc) = mean( VAFc ); 
   VAFmin(appl,lc)  = min( VAFc ); 
   %
   clear data K model sys u y ye VAFc
end
%
format short e
disp( ' ' );
disp( 'Execution time (seconds)' )
disp( ' ' )
%
% Timing
%
disp( '-------------------------------------------' )
disp( '  App.    Dimensions  Execution time (sec.)' )
disp( '-------------------------------------------' )
disp( '  #     t  m  l  s  n   slmoen4   n4sid    ' )
disp( '-------------------------------------------' )
for L = 1 : appl
   disp( [sprintf('%3d',dim(L,1)),  sprintf('%6d',dim(L,2)),...
          sprintf('%3d',dim(L,3:6)),sprintf('  %7.2f',timing(L,:))] )
end   
disp( '-------------------------------------------' )
%
disp( ' ' )
disp( 'Relative errors, norm(y - ye,1)/norm(y,1)' )
disp( ' ' )
%
% Errors
%
disp( '-------------------------------------------' )
disp( '  App.    Dimensions     Relative errors   ' )
disp( '-------------------------------------------' )
disp( '  #     t  m  l  s  n   slmoen4     n4sid  ' )
disp( '-------------------------------------------' )
for L = 1 : appl
   disp( [sprintf('%3d',dim(L,1)),  sprintf('%6d',dim(L,2)),...
          sprintf('%3d',dim(L,3:6)),sprintf('  %0.2e',err(L,:))] )
end   
disp( '-------------------------------------------' )
%
disp( ' ' )
disp( 'Means of Variance-Accounted-For (in percentages), mean(VAF),' )
disp( 'for all outputs of a system ' )
disp( ' ' )
%
% Mean(VAF)
%
disp( '-------------------------------------------' )
disp( '  App.    Dimensions        mean(VAF)      ' )
disp( '-------------------------------------------' )
disp( '  #     t  m  l  s  n   slmoen4   n4sid    ' )
disp( '-------------------------------------------' )
for L = 1 : appl
   disp( [sprintf('%3d',dim(L,1)),  sprintf('%6d',dim(L,2)),...
          sprintf('%3d',dim(L,3:6)),sprintf('  %7.2f',VAFmean(L,:))] )
end   
disp( '-------------------------------------------' )
%
disp( ' ' )
disp( 'Mimimum of Variance-Accounted-For (in percentages), min(VAF),' )
disp( 'for all outputs of a system ' )
disp( ' ' )
%
% Min(VAF)
%
disp( '-------------------------------------------' )
disp( '  App.    Dimensions        min(VAF)       ' )
disp( '-------------------------------------------' )
disp( '  #     t  m  l  s  n   slmoen4   n4sid    ' )
disp( '-------------------------------------------' )
for L = 1 : appl
   disp( [sprintf('%3d',dim(L,1)),  sprintf('%6d',dim(L,2)),...
          sprintf('%3d',dim(L,3:6)),sprintf('  %7.2f',VAFmin(L,:))] )
end   
disp( '-------------------------------------------' )
%
% Graphical representations.
%
disp( ' ' )
disp( 'Graphical representations' )
disp( ' ' )
disp( 'SLICOT function slmoen4, algorithm alg, versus n4sid with options' )
disp( '''N4Weight'' = ''MOESP'', and ''order'', ''CovarianceMatrix'', and ''N4Horizon'' = ''Auto'',' )
disp( 'set as required by given_n, covr, and given_N4H, respectively.' )
disp( ' ' )
disp( 'Timings' )
disp( ' ' )
%
timesC = timing;
%
set( axes, 'FontSize', 16 )
axis on 
bar( timesC )
legend( 'slmoen4', 'n4sid', 'location', 'Best' )
xlabel( 'Application #' )
ylabel( 'Time (sec.)' )
disp( 'Press any key to continue' )
disp( ' ' )
pause
%
% Speed-up SLICOT slmoen4, algorithm alg, over n4sid with option 'MOESP',
% and options for 'order', 'CovarianceMatrix', and 'N4Horizon' set as
% required by given_n, covr, and given_N4H, respectively.
%
disp( 'Speed-up factor' )
disp( ' ' )
timtmp = timing(:,1);
timtmp( timtmp == 0) = .01;
%
tim2M = timing(:,2)./timtmp;
%
figure
set( axes, 'FontSize', 16 )
axis on 
bar( tim2M )
xlabel( 'Application #' )
ylabel( 'Speed-up factor' )
disp( 'Press any key to continue' )
disp( ' ' )
pause
%
% Relative errors.
%
disp( 'Relative errors' )
disp( ' ' )
figure
set( axes, 'FontSize', 16 )
axis on 
bar( err )
legend( 'slmoen4', 'n4sid', 'location', 'Best' )
xlabel( 'Application #' )
ylabel( 'Relative errors' )
disp( 'Press any key to continue' )
disp( ' ' )
pause
%
% end of cmp_ident.m
