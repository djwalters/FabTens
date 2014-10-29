function [EStarRaw, GStarRaw, nuStarRaw, EStarAr, GStarAr, nuStarAr,...
    YoungsRaw, LYoungsRaw, UYoungsRaw,...
    YoungsAr, LYoungsAr, UYoungsAr,...
    YoungsArTr, LYoungsArTr, UYoungsArTr,...
    ShearRaw, LShearRaw, UShearRaw,...
    ShearAr, LShearAr, UShearAr,...
    ShearArTr, LShearArTr, UShearArTr,...
    PoissonRaw, LPoissonRaw, UPoissonRaw,...
    PoissonAr, LPoissonAr, UPoissonAr,...
    PoissonArTr, LPoissonArTr, UPoissonArTr,...
    SRaw, SAr, SArTr] = ...
    MechModuli(phii, N3, rho_hat, R_hat, C, Unc,Ar,Tr)
% MechModuli.m
% function [EStarRaw, GStarRaw, nuStarRaw, EStarAr, GStarAr, nuStarAr,...
%     YoungsRaw, LYoungsRaw, UYoungsRaw,...
%     YoungsAr, LYoungsAr, UYoungsAr,...
%     YoungsArTr, LYoungsArTr, UYoungsArTr,...
%     ShearRaw, LShearRaw, UShearRaw,...
%     ShearAr, LShearAr, UShearAr,...
%     ShearArTr, LShearArTr, UShearArTr,...
%     PoissonRaw, LPoissonRaw, UPoissonRaw,...
%     PoissonAr, LPoissonAr, UPoissonAr,...
%     PoissonArTr, LPoissonArTr, UPoissonArTr,...
%     SRaw] = ...
%     MechModuli(phii, N3, rho_hat, R_hat, C, Unc,Ar,Tr)
% This funtion calculates the effective Young's moduli and Shear Moduli of
% the sample.  These are calculated using stereologically determined
% isotropic constants (E* and G*) and the calculated contact tensor
% coefficients.
%
% INPUTS:
%   phii - Solid volume fraction from SkyScan CTAn.
%
%   N3 - 3-D Coordination number calculated from Edens' SSAnalyzer or
%   appropriate parameter from SkyScan CTAn software.
%
%   rho_hat - 3-D bond radius
%
%   R_hat - 3-D grain radius
%
%   C - Contact tensor coefficients.
%           Format: C(row,column,run number)
%           i.e. k = 3, run number = 3
%                   C = [F11,F12,F13;
%                        F21,F22,F23;
%                        F31,F32,F33]
%                   F22 = row 2, columnn 2, run 3 = C(2,2,3)
%
%   Unc - Confidence intervals for the contact tensor coefficients.
%           Format: Unc(row, column, run number, lower(1)/upper(2)
%           interval)
%
%   Ar - Grain shape factor aspect ratio
%
%   Tr - Tensor ratio (between MIL Tensor and Contact Tensor)
%
% OUTPUTS:
%   EStar(---) - Scalar Young's modulus calculated from measurable
%   microstructural parameters.
%       *Raw - Results produced by the original contact tensor model.
%       *Ar - Results as produced by incorporating only the scalar shape
%       factor aspect ratio, Ar.
%
%   GStar(---) - Scalar shear modulus calculated from measurable
%   microstructural parameters.
%
%   nuStar(---) - Scalar Poisson's ratio calculated from measurable
%   microstructural parameters.
%
%   Youngs(---) - A structure of anisotropic Young's moduli along each of the 
%   principal directions of contact (i.e [E1,E2,E3]).
%       *Raw - Results as produced by the orginal contact tensor model.
%       *Ar - Results as produced by incorporating only the scalar shape
%       factor aspect ratio, Ar.
%       *ArTr - Results as produced by incorporating both the scalar shape
%       factor aspect ratio, Ar, and the tensor ratio, Tr.
%   
%   LYoungs(---)/UYoungs(---) - Lower and upper confidence intervals for
%   Youngs---.
%
%   Shear(---) - A structure of anisotropic shear moduli along aligning with 
%   the principal directions of contact (i.e. [G23,G13,G12]).
%
%   LShear(---)/UShear(---) - Lower and upper confidence intervals for
%   Shear(---).
%
%   Poisson(---) - A structure of anisotropic Poisson's ratio alignining with 
%   the principal directions of contact (i.e. [nu23,nu13,nu12]).
%
%   LPossion(---)/UPoisson(---) - Lower and upper confidence intervals for
%   Poisson(---).
%
%   S(---) - A 6x6 array that represents the anisotropic compliance matrix.
%   Its maximum deviation from isotropy is orthotropic symmetry (9 
%   independent coefficients) as restricted by the use of a 2nd order
%   contact tensor.
%
%   v1.0, June 17, 2013
%   v1.1, June 21, 2013, Added data for plotting moduli and confidence
%   intervals
%   v1.2, October 29, 2014, Added shape factor correction and appropriate
%   outputs for gains made by shape factors
% AUTHOR: David J. Walters; Montana State University

%% Initialize variables internal to function
% E_ice = 9330e6;     %Young's modulus ice (Petrenko and Whitworth 1999)
% G_ice = 3520e6;     %Shear modulus ice (Petrenko and Whitworth 1999)
E_ice = 1000e6;
G_ice = 377e6;
xi = G_ice/E_ice;   %Stiffness ratio
f = 1/3;            %Isotropic value of contact tensor

%% Calculate scalar (isotropic) moduli from stereological data
EStarRaw = (1/2) * (rho_hat/R_hat)^2 * phii * N3 * E_ice * (2+3*xi)/(4+xi);
GStarRaw = (1/40) * (rho_hat/R_hat)^2 * phii * N3 * E_ice * (2+3*xi);
nuStarRaw = (1-xi)/(4+xi);

EStarAr = EStarRaw * Ar;
GStarAr = GStarRaw * Ar;
nuStarAr = nuStarRaw;

%% Calculate anisotropic moduli as modified by the contact tensor

SRaw = zeros(6,6);
SAr = zeros(6,6);
SArTr = zeros(6,6);
LSRaw = zeros(6,6);
LSAr = zeros(6,6);
LSArTr = zeros(6,6);
USRaw = zeros(6,6);
USAr = zeros(6,6);
USArTr = zeros(6,6);
ERaw = zeros(3);
EAr = zeros(3);
EArTr = zeros(3);
LERaw = zeros(3);
LEAr = zeros(3);
LEArTr = zeros(3);
UERaw = zeros(3);
UEAr = zeros(3);
UEArTr = zeros(3);
nuRaw = zeros(3);
nuAr = zeros(3);
nuArTr = zeros(3);
LnuRaw = zeros(3);
LnuAr = zeros(3);
LnuArTr = zeros(3);
UnuRaw = zeros(3);
UnuAr = zeros(3);
UnuArTr = zeros(3);

% Work on upper left 3x3 corner of compliance matrix
k = 0;
for i = 1:3
    for j = 1:3
        if i == j % Diagonal values
            SRaw(i,j) = (1/EStarRaw) * (f/C(i,j,:))^2;
            SAr(i,j) = (1/EStarAr) * (f/C(i,j,:))^2;
            SArTr(i,j) = (1/EStarAr) * (f/C(i,j,:))^2 * (1/Tr(i))^2;
            
            LSRaw(i,j) = (1/EStarRaw) * (f/Unc(i,j,:,1))^2;  % Lower 95% CI
            LSAr(i,j) = (1/EStarAr) * (f/Unc(i,j,:,1))^2;
            LSArTr(i,j) = (1/EStarAr) * (f/Unc(i,j,:,1))^2 * (1/Tr(i))^2;
            
            USRaw(i,j) = (1/EStarRaw) * (f/Unc(i,j,:,2))^2;  % Upper 95% CI
            USAr(i,j) = (1/EStarAr) * (f/Unc(i,j,:,2))^2;
            USArTr(i,j) = (1/EStarAr) * (f/Unc(i,j,:,2))^2 * (1/Tr(i))^2;
            
            ERaw(i) = SRaw(i,j)^(-1);
            EAr(i) = SAr(i,j)^(-1);
            EArTr(i) = SArTr(i,j)^(-1);
            
            LERaw(i) = LSRaw(i,j)^(-1);   % Lower 95% CI
            LEAr(i) = LSAr(i,j)^(-1);
            LEArTr(i) = LSArTr(i,j)^(-1);
            
            UERaw(i) = USRaw(i,j)^(-1);   % Upper 95% CI
            UEAr(i) = USAr(i,j)^(-1);
            UEArTr(i) = USArTr(i,j)^(-1);
        else
            SRaw(i,j) = (nuStarRaw/EStarRaw) * (f^2/(C(j,j,:)*C(i,i,:)));
            SAr(i,j) = (nuStarAr/EStarAr) * (f^2/(C(j,j,:)*C(i,i,:)));
            SArTr(i,j) = (nuStarAr/EStarAr) * (f^2/(C(j,j,:)*C(i,i,:))) * (1/(Tr(i)*Tr(j)));
            
            LSRaw(i,j) = (nuStarRaw/EStarRaw) * (f^2/(Unc(j,j,:,1)*Unc(i,i,:,1)));
            LSAr(i,j) = (nuStarAr/EStarAr) * (f^2/(Unc(j,j,:,1)*Unc(i,i,:,1)));
            LSArTr(i,j) = (nuStarAr/EStarAr) * (f^2/(Unc(j,j,:,1)*Unc(i,i,:,1))) * (1/(Tr(i)*Tr(j)));
            
            USRaw(i,j) = (nuStarRaw/EStarRaw) * (f^2/(Unc(j,j,:,2)*Unc(i,i,:,2)));
            USAr(i,j) = (nuStarAr/EStarAr) * (f^2/(Unc(j,j,:,2)*Unc(i,i,:,2)));
            USArTr(i,j) = (nuStarAr/EStarAr) * (f^2/(Unc(j,j,:,2)*Unc(i,i,:,2))) * (1/(Tr(i)*Tr(j)));
            
            if j > i
                k = k + 1;
            nuRaw(k) = SRaw(i,j)*ERaw(i);
            nuAr(k) = SAr(i,j)*EAr(i);
            nuArTr(k) = SArTr(i,j)*EArTr(i);
            
            LnuRaw(k) = LSRaw(i,j)*ERaw(i);
            LnuAr(k) = LSAr(i,j)*EAr(i);
            LnuArTr(k) = LSArTr(i,j)*EArTr(i);
            
            UnuRaw(k) = USRaw(i,j)*ERaw(i);
            UnuAr(k) = USAr(i,j)*EAr(i);
            UnuArTr(k) = USArTr(i,j)*EArTr(i);
            end
        end
    end
end

% Shear values
SRaw(4,4) = (1/(2*GStarRaw)) * (f^2/(C(2,2,:)*C(3,3,:)));
    LSRaw(4,4) = (1/(2*GStarRaw)) * (f^2/(Unc(2,2,:,1)*Unc(3,3,:,1)));
    USRaw(4,4) = (1/(2*GStarRaw)) * (f^2/(Unc(2,2,:,2)*Unc(3,3,:,2)));
SAr(4,4) = (1/(2*GStarAr)) * (f^2/(C(2,2,:)*C(3,3,:)));
    LSAr(4,4) = (1/(2*GStarAr)) * (f^2/(Unc(2,2,:,1)*Unc(3,3,:,1)));
    USAr(4,4) = (1/(2*GStarAr)) * (f^2/(Unc(2,2,:,2)*Unc(3,3,:,2)));
SArTr(4,4) = (1/(2*GStarAr)) * (f^2/(C(2,2,:)*C(3,3,:))) * (1/(Tr(2)*Tr(3)));
    LSArTr(4,4) = (1/(2*GStarAr)) * (f^2/(Unc(2,2,:,1)*Unc(3,3,:,1))) * (1/(Tr(2)*Tr(3)));
    USArTr(4,4) = (1/(2*GStarAr)) * (f^2/(Unc(2,2,:,2)*Unc(3,3,:,2))) * (1/(Tr(2)*Tr(3)));
    
SRaw(5,5) = (1/(2*GStarRaw)) * (f^2/(C(1,1,:)*C(3,3,:)));
    LSRaw(5,5) = (1/(2*GStarRaw)) * (f^2/(Unc(1,1,:,1)*Unc(3,3,:,1)));
    USRaw(5,5) = (1/(2*GStarRaw)) * (f^2/(Unc(1,1,:,2)*Unc(3,3,:,2)));
SAr(5,5) = (1/(2*GStarAr)) * (f^2/(C(1,1,:)*C(3,3,:)));
    LSAr(5,5) = (1/(2*GStarAr)) * (f^2/(Unc(1,1,:,1)*Unc(3,3,:,1)));
    USAr(5,5) = (1/(2*GStarAr)) * (f^2/(Unc(1,1,:,2)*Unc(3,3,:,2)));
SArTr(5,5) = (1/(2*GStarAr)) * (f^2/(C(1,1,:)*C(3,3,:))) * (1/(Tr(1)*Tr(3)));
    LSArTr(5,5) = (1/(2*GStarAr)) * (f^2/(Unc(1,1,:,1)*Unc(3,3,:,1))) * (1/(Tr(1)*Tr(3)));
    USArTr(5,5) = (1/(2*GStarAr)) * (f^2/(Unc(1,1,:,2)*Unc(3,3,:,2))) * (1/(Tr(1)*Tr(3)));
    
SRaw(6,6) = (1/(2*GStarRaw)) * (f^2/(C(1,1,:)*C(2,2,:)));
    LSRaw(6,6) = (1/(2*GStarRaw)) * (f^2/(Unc(1,1,:,1)*Unc(2,2,:,1)));
    USRaw(6,6) = (1/(2*GStarRaw)) * (f^2/(Unc(1,1,:,2)*Unc(2,2,:,2)));
SAr(6,6) = (1/(2*GStarAr)) * (f^2/(C(1,1,:)*C(2,2,:)));
    LSAr(6,6) = (1/(2*GStarAr)) * (f^2/(Unc(1,1,:,1)*Unc(2,2,:,1)));
    USAr(6,6) = (1/(2*GStarAr)) * (f^2/(Unc(1,1,:,2)*Unc(2,2,:,2)));
SArTr(6,6) = (1/(2*GStarAr)) * (f^2/(C(1,1,:)*C(2,2,:))) * (1/(Tr(1)*Tr(2)));
    LSArTr(6,6) = (1/(2*GStarAr)) * (f^2/(Unc(1,1,:,1)*Unc(2,2,:,1))) * (1/(Tr(1)*Tr(2)));
    USArTr(6,6) = (1/(2*GStarAr)) * (f^2/(Unc(1,1,:,2)*Unc(2,2,:,2))) * (1/(Tr(1)*Tr(2)));
    
% Raw Outputs
YoungsRaw.x = ERaw(1);
    LYoungsRaw.x = LERaw(1);
    UYoungsRaw.x = UERaw(1);
YoungsRaw.y = ERaw(2);
    LYoungsRaw.y = LERaw(2);
    UYoungsRaw.y = UERaw(2);
YoungsRaw.z = ERaw(3);
    LYoungsRaw.z = LERaw(3);
    UYoungsRaw.z = UERaw(3);

ShearRaw.xy = SRaw(6,6)^(-1)/2;
    LShearRaw.xy = LSRaw(6,6)^(-1);
    UShearRaw.xy = USRaw(6,6)^(-1);
ShearRaw.xz = SRaw(5,5)^(-1)/2;
    LShearRaw.xz = LSRaw(5,5)^(-1);
    UShearRaw.xz = USRaw(5,5)^(-1);
ShearRaw.yz = SRaw(4,4)^(-1)/2;
    LShearRaw.yz = LSRaw(4,4)^(-1);
    UShearRaw.yz = USRaw(4,4)^(-1);

PoissonRaw.xy = nuRaw(1);
    LPoissonRaw.xy = LnuRaw(1);
    UPoissonRaw.xy = UnuRaw(1);
PoissonRaw.xz = nuRaw(2);
    LPoissonRaw.xz = LnuRaw(2);
    UPoissonRaw.xz = UnuRaw(2);
PoissonRaw.yz = nuRaw(3);
    LPoissonRaw.yz = LnuRaw(3);
    UPoissonRaw.yz = UnuRaw(3);

% Ar adjusted outputs
YoungsAr.x = EAr(1);
    LYoungsAr.x = LEAr(1);
    UYoungsAr.x = UEAr(1);
YoungsAr.y = EAr(2);
    LYoungsAr.y = LEAr(2);
    UYoungsAr.y = UEAr(2);
YoungsAr.z = EAr(3);
    LYoungsAr.z = LEAr(3);
    UYoungsAr.z = UEAr(3);

ShearAr.xy = SAr(6,6)^(-1)/2;
    LShearAr.xy = LSAr(6,6)^(-1);
    UShearAr.xy = USAr(6,6)^(-1);
ShearAr.xz = SAr(5,5)^(-1)/2;
    LShearAr.xz = LSAr(5,5)^(-1);
    UShearAr.xz = USAr(5,5)^(-1);
ShearAr.yz = SAr(4,4)^(-1)/2;
    LShearAr.yz = LSAr(4,4)^(-1);
    UShearAr.yz = USAr(4,4)^(-1);

PoissonAr.xy = nuAr(1);
    LPoissonAr.xy = LnuAr(1);
    UPoissonAr.xy = UnuAr(1);
PoissonAr.xz = nuAr(2);
    LPoissonAr.xz = LnuAr(2);
    UPoissonAr.xz = UnuAr(2);
PoissonAr.yz = nuAr(3);
    LPoissonAr.yz = LnuAr(3);
    UPoissonAr.yz = UnuAr(3);
 
% Ar and Tr adjusted outputs
YoungsArTr.x = EArTr(1);
    LYoungsArTr.x = LEArTr(1);
    UYoungsArTr.x = UEArTr(1);
YoungsArTr.y = EArTr(2);
    LYoungsArTr.y = LEArTr(2);
    UYoungsArTr.y = UEArTr(2);
YoungsArTr.z = EArTr(3);
    LYoungsArTr.z = LEArTr(3);
    UYoungsArTr.z = UEArTr(3);

ShearArTr.xy = SArTr(6,6)^(-1)/2;
    LShearArTr.xy = LSArTr(6,6)^(-1);
    UShearArTr.xy = USArTr(6,6)^(-1);
ShearArTr.xz = SArTr(5,5)^(-1)/2;
    LShearArTr.xz = LSArTr(5,5)^(-1);
    UShearArTr.xz = USArTr(5,5)^(-1);
ShearArTr.yz = SArTr(4,4)^(-1)/2;
    LShearArTr.yz = LSArTr(4,4)^(-1);
    UShearArTr.yz = USArTr(4,4)^(-1);

PoissonArTr.xy = nuArTr(1);
    LPoissonArTr.xy = LnuArTr(1);
    UPoissonArTr.xy = UnuArTr(1);
PoissonArTr.xz = nuArTr(2);
    LPoissonArTr.xz = LnuArTr(2);
    UPoissonArTr.xz = UnuArTr(2);
PoissonArTr.yz = nuArTr(3);
    LPoissonArTr.yz = LnuArTr(3);
    UPoissonArTr.yz = UnuArTr(3);

end