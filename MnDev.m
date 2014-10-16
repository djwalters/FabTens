function [Dev2,Dev4] = MnDev(F2,F4)
% This function computes the mean and deviatoric parts of the symmetric 2nd
% and 4th rank contact tensors (arrays)
%
%Input:
%   F2 - symmetric, 2nd order contact tensor
%   F4 - symmetric, 4th order contact tensor
%   fn - field names, "xy", "xz", "yz"
%   n - number of fields
%
%Output:
%   Mn - the scalar, mean value of input N2
%   Dev2 - the deviatoric part of N2, a 2nd rank tensor
%   Dev4 - the deviatoric part of N4, a 4th rank tensor

%% Intialize variables
dim = 3;
I = eye(dim);   % Define 2nd order identity/kronecker tensor

%% Calculate mean and deviatoric parts

% Calculate 2nd order mean
mn = trace(F2) / dim;
Mn = trace(F2);     %Density Distribution Mean is 1st invariant of N2

% Calculate 2nd order deviator
DevN = F2 - mn*I;

% Calculate 4th order deviator
for i = 1:dim
    
    for j = 1:dim
        for k = 1:dim
            for l = 1:dim
                K(i,j,k,l) = 1/3 * (I(i,j) * I(k,l) +...
                    I(i,k) * I(j,l )+...
                    I(i,l) * I(j,k));
                A(i,j,k,l) = 1/6 * (I(i,j) * F2(k,l) +...
                    I(k,l) * F2(i,j) +...
                    I(i,k) * F2(j,l) +...
                    I(i,l) * F2(j,k) +...
                    I(j,k) * F2(i,l) +...
                    I(j,l) * F2(i,k));
            end
        end
    end
end

% See Equations 2.20 and 2.21 in Shertzer's dissertation for
% explanation on the following.  See Kanatani (1984)
Dev2 = 15/2 * DevN;
Dev4= 315/8 * (F4 - 6/7*A + 3/35*K);


end

%EOF