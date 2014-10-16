function [E_star, G_star, nu_star, Youngs, LYoungs, UYoungs,...
    Shear, LShear, UShear, Poisson, LPoisson, UPoisson, S] = ...
    MechModuli(phii, N3, rho_hat, R_hat, C, Unc)
% This funtion calculates the effective Young's moduli and Shear Moduli of
% the sample.  These are calculated using stereologically determined
% isotropic constants (E* and G*) and the calculated contact tensor
% coefficients.
%
% INPUTS:
%   phii - Solid volume fraction (currently set to the xy plane
%               calculated from Edens' SSAnalyzer -- Needs to be replaced
%               with 3-D value from SkyScan CTAn software).
%   N3 - 3-D Coordination number calculated from Edens' SSAnalyzer or
%           appropriate parameter from SkyScan CTAn software.
%   rho_hat - 3-D bond radius
%   R_hat - 3-D grain radius
%   C - Contact tensor coefficients.
%           Format: C(row,column,run number)
%           i.e. k = 3, run number = 3
%                   C = [F11,F12,F13;
%                        F21,F22,F23;
%                        F31,F32,F33]
%                   F22 = row 2, columnn 2, run 3 = C(2,2,3)
%
% OUTPUTS:
%   Youngs - A structure of anisotropic Young's moduli along each of the 
%               principal directions of contact (i.e [E1,E2,E3]).
%   Shear - A structure of anisotropic shear moduli along aligning with the 
%              principal directions of contact (i.e. [G23,G13,G12]).
%   Poisson - A structure of anisotropic Poisson's ratio alignining with the
%              principal directions of contact (i.e. [nu23,nu13,nu12]).
%   S - A 6x6 array that represents the anisotropic compliance matrix.
%          Its maximum deviation from isotropy is orthotropic symmetry (9 
%          independent coefficients) as restricted by the use of a 2nd
%          order contact tensor.
%
% Function Author: 
%   David J. Walters
%   Montana State University
%   walters.david.j@gmail.com, david.walters@ce.montana.edu
%   406-924-9597
%   v1.0, June 17, 2013
%   v1.1, June 21, 2013, Added data for plotting moduli and confidence
%   intervals

%% Initialize variables internal to function
E_ice = 9330e6;     %Young's modulus ice (Petrenko and Whitworth 1999)
G_ice = 3520e6;     %Shear modulus ice (Petrenko and Whitworth 1999)
xi = G_ice/E_ice;   %Stiffness ratio
f = 1/3;            %Isotropic value of contact tensor

%% Calculate scalar (isotropic) moduli from stereological data
E_star = (1/2) * (rho_hat/R_hat)^2 * phii * N3 * E_ice * (2+3*xi)/(4+xi);
G_star = (1/40) * (rho_hat/R_hat)^2 * phii * N3 * E_ice * (2+3*xi);
nu_star = (1-xi)/(4+xi);

%% Calculate anisotropic moduli as modified by the contact tensor

S = zeros(6,6);
LS = zeros(6,6);
US = zeros(6,6);
E = zeros(3);
LE = zeros(3);
UE = zeros(3);
nu = zeros(3);
Lnu = zeros(3);
Unu = zeros(3);

% Work on upper left 3x3 corner of compliance matrix
k = 0;
for i = 1:3
    for j = 1:3
        if i == j % Diagonal values
            S(i,j) = (1/E_star) * (f/C(i,j,:))^2;
            LS(i,j) = (1/E_star) * (f/Unc(i,j,:,1))^2;  % Lower 95% CI
            US(i,j) = (1/E_star) * (f/Unc(i,j,:,2))^2;  % Upper 95% CI
            E(i) = S(i,j)^(-1);
            LE(i) = LS(i,j)^(-1);   % Lower 95% CI
            UE(i) = US(i,j)^(-1);   % Upper 95% CI
        else
            S(i,j) = (nu_star/E_star) * (f^2/(C(j,j,:)*C(i,i,:)));
            LS(i,j) = (nu_star/E_star) * (f^2/(Unc(j,j,:,1)*Unc(i,i,:,1)));
            US(i,j) = (nu_star/E_star) * (f^2/(Unc(j,j,:,2)*Unc(i,i,:,2)));
            if j > i
                k = k + 1;
            nu(k) = S(i,j)*E(i);
            Lnu(k) = LS(i,j)*E(i);
            Unu(k) = US(i,j)*E(i);
            end
        end
    end
end

% Shear values
S(4,4) = (1/(2*G_star)) * (f^2/(C(2,2,:)*C(3,3,:)));
    LS(4,4) = (1/(2*G_star)) * (f^2/(Unc(2,2,:,1)*Unc(3,3,:,1)));
    US(4,4) = (1/(2*G_star)) * (f^2/(Unc(2,2,:,2)*Unc(3,3,:,2)));
S(5,5) = (1/(2*G_star)) * (f^2/(C(1,1,:)*C(3,3,:)));
    LS(5,5) = (1/(2*G_star)) * (f^2/(Unc(1,1,:,1)*Unc(3,3,:,1)));
    US(5,5) = (1/(2*G_star)) * (f^2/(Unc(1,1,:,2)*Unc(3,3,:,2)));
S(6,6) = (1/(2*G_star)) * (f^2/(C(1,1,:)*C(2,2,:)));
    LS(6,6) = (1/(2*G_star)) * (f^2/(Unc(1,1,:,1)*Unc(2,2,:,1)));
    US(6,6) = (1/(2*G_star)) * (f^2/(Unc(1,1,:,2)*Unc(2,2,:,2)));

% Outputs
Youngs.x = E(1);
    LYoungs.x = LE(1);
    UYoungs.x = UE(1);
Youngs.y = E(2);
    LYoungs.y = LE(2);
    UYoungs.y = UE(2);
Youngs.z = E(3);
    LYoungs.z = LE(3);
    UYoungs.z = UE(3);


Shear.xy = S(6,6)^(-1);
    LShear.xy = LS(6,6)^(-1);
    UShear.xy = US(6,6)^(-1);
Shear.xz = S(5,5)^(-1);
    LShear.xz = LS(5,5)^(-1);
    UShear.xz = US(5,5)^(-1);
Shear.yz = S(4,4)^(-1);
    LShear.yz = LS(4,4)^(-1);
    UShear.yz = US(4,4)^(-1);

Poisson.xy = nu(1);
    LPoisson.xy = Lnu(1);
    UPoisson.xy = Unu(1);
Poisson.xz = nu(2);
    LPoisson.xz = Lnu(2);
    UPoisson.xz = Unu(2);
Poisson.yz = nu(3);
    LPoisson.yz = Lnu(3);
    UPoisson.yz = Unu(3);

end