function v = vaf(y,y_est)
% vaf        Compute the percentage Variance Accounted For (VAF)
%            between two signals. the VAF is calculated as:
%                        variance(y-y_est)
%               v= (1 -  ----------------  ) * 100%
%           	             variance(y)
% 
%            The VAF of two signals that are the same is
%            100%. If they differ, the  VAF will be lower. 
%            When y and y_est have multiple columns, the VAF 
%            is calculated for every column in y and y_est. 
%            The VAF is often used to verify the
%            correctness of a model, by comparing the real
%            output with the estimated output of the model.
% 
% Syntax:
%            v = vaf(y,y_estimate)
% Input:     
%  y         Signal 1, often the real output.
%  y_est     Signal 2, often the estimated output of a model.
%
% Output:
%   v        VAF, computed for the two signals

% RELEASE 2.0 of SLICOT System Identification Toolbox.
%
% Bert Haverkamp, April 1996
% copyright (c) 1996 B.R.J. Haverkamp

% Revised Vasile Sima, December 1999

if nargin<2
  error('Not enough input variables')
end
if size(y,2)>size(y,1)
  y=y';
end
if size(y_est,2)>size(y_est,1)
  y_est=y_est';
end

N=size(y,1);

if ~(size(y_est,1)==N)
  error('Both signals should have same length')
end

v=diag(100*(eye(size(y,2))-cov(y-y_est)./cov(y)));




