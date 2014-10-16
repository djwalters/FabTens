function [zcrit,z2,z4] = ConvgTest(N,dim,D2,D4,alpha,fn,n)
% This function performs a hypothesis test to determine whether the higher
% rank tensors are required to more accurately describe the distribution
% density
%
% Input:
%   N - # of samples (contact normal vectors)
%   dim - dimension of individual tensors; 2 for 2-D problem, 3 for 3-D 
%       problem
%   D2 - symmetric, 2nd rank deviatoric tensor
%   D4 - symmetric, 4th rank deviatoric tensor
%   alpha - user defined significance of hypothesis test
%   fn - field names of the three orthogonal planes: "xy", "xz", "yz"
%   n - dimension of problem, e.g. n = 3 for 3-D
%
% Output:
%   zcrit - critical value
%   z2 - test statistic for 2nd rank tensor
%   z4 - test statistic for 4th rank tensor
%
% Information
%   1.) Null hypothesis: true distribution has anisotropy characterized by
%           D2 or D4.
%
%   2.) Alternative hypothesis: true distribution is described by a lesser
%           order tensor -- D2 not required, scalar works, or D4 not
%           required, use D2
%
%   3.) Calculate test statistic: liklihood ratio --  Probability that data
%           comes from uniform distribution -or- probability that data
%           comes from a "D2" distribution


z2 = N/(dim^2) * sum(sum(D2.^2));
z4 = [];     %for now...

% Lookup critical value: liklihood ratio behaves according to
% Chi-square
zcrit = chi2inv(1-alpha,2);

% Compare test statistic to critical value and repeat for 4th rank if
% necessary
if z2 <= zcrit
    disp('D2 and D4 not required')
elseif z2 > zcrit
    disp(['** D2 Required ** - '])
    z4 = N/(dim^4) * sum(sum(sum(sum(D4.^2))));
    if z4 <= zcrit
        disp('D4 not required')
    elseif z4 > zcrit
        disp(['** D4 Required ** - '])
    end
end

end

%EOF