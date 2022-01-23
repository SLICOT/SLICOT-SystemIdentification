function y = NNout( n, l, wb, x )
%NNout  Computes the output of a set of neural networks with the structure
%
%             - tanh(w1'*x+b1) -
%           /      :             \
%         x ---    :           --- sum(ws(i)*...)+ b(n+1)  --- y,
%           \      :             /
%             - tanh(wn'*x+bn) -
%
%       where x, w1, ..., wn are vectors of length m, ws is a vector
%       of length n, b(1), ..., b(n+1) are scalars, and n is the
%       number of neurons in the hidden layer.
%       Such a network is used for each output variables.
%
%       Y = NNOUT(N,L,WB,X)  computes the output Y of a set of L neural
%       networks with N neurons, given the input trajectory X, 
%       X = [x1'; x2'; ..., xT']', with each xj a vector of length M,
%       and the parameter vector WB, containing the weights and biases 
%       of the network. The vector WB is partitioned into L vectors of
%       length N*(M+2)+1, WB = [ wb(1), ..., wb(L) ]. Each wb(k),  k = 1, 
%       ..., L, corresponds to one output variable, and is defined as
%       wb(k) = [ w1(1), ..., w1(M), ..., wn(1), ..., wn(M),
%                 ws(1), ..., ws(N), b(1), ..., b(N+1) ],
%       where wi(j) are the weights of the hidden layer,
%       ws(i) are the weights of the linear output layer, and
%       b(i) are the biases, as in the scheme above.

%   RELEASE 2.0 of SLICOT System Identification Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001.
%
%   Revisions:
%   -
%

% Check the length of input data.

if nargin < 4,
   error( 'Wrong number of input arguments (should be at least 4)' )
elseif nargout < 1,
   error( 'Wrong number of output arguments (should be at least 1)' )
end

[ t, m ] = size( x );
lwb = l*( n*( m + 2 ) + 1 );
if length( wb ) < lwb,
    error( [ 'WB is not large enough (it should be at least ', num2str( lwb ),')' ] )
end
if l > 0,  lwb = lwb/l;  end  
y = zeros(t,l);

lk = 0;  is = n*m;  ib = n*( m + 1 );
for k = 0 : l - 1,
    bnp1 = wb(ib+n+1+lk);
    W = wb( 1+lk : n*m+lk );  W = reshape( W, m, n );
    for i = 1 : t,
        if n > 0,
           w = W'*x(i,:)' + wb( ib+1+lk : ib+n+lk );
        else
           w = W'*x(i,:)';
        end           
        w = tanh(w);
        if n > 0,
           y(i,k+1) = w'*wb( is+1+lk : is+n+lk ) + bnp1;
        else
           y(i,k+1) = bnp1;
        end           
    end
    lk = lk + lwb;
end
%
% end NNout

