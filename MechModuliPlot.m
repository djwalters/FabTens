function MechModuliPlot( EStarRaw, GStarRaw, nuStarRaw, EStarAr, GStarAr, nuStarAr,...
    YoungsRaw, LYoungsRaw, UYoungsRaw,...
    YoungsAr, LYoungsAr, UYoungsAr,...
    YoungsArTr, LYoungsArTr, UYoungsArTr,...
    ShearRaw, LShearRaw, UShearRaw,...
    ShearAr, LShearAr, UShearAr,...
    ShearArTr, LShearArTr, UShearArTr,...
    PoissonRaw, LPoissonRaw, UPoissonRaw,...
    PoissonAr, LPoissonAr, UPoissonAr,...
    PoissonArTr, LPoissonArTr, UPoissonArTr,...
    SRaw, SAr, SArTr, endtime, idx)
% MechModuliPlot.m
% function MechModuliPlot( EStarRaw, GStarRaw, nuStarRaw, EStarAr, GStarAr, nuStarAr,...
%     YoungsRaw, LYoungsRaw, UYoungsRaw,...
%     YoungsAr, LYoungsAr, UYoungsAr,...
%     YoungsArTr, LYoungsArTr, UYoungsArTr,...
%     ShearRaw, LShearRaw, UShearRaw,...
%     ShearAr, LShearAr, UShearAr,...
%     ShearArTr, LShearArTr, UShearArTr,...
%     PoissonRaw, LPoissonRaw, UPoissonRaw,...
%     PoissonAr, LPoissonAr, UPoissonAr,...
%     PoissonArTr, LPoissonArTr, UPoissonArTr,...
%     SRaw, SAr, SArTr )
%
% This function plots the results of the function MechModuli.m.  It should
% be able to appropriately handle time series data.
%
% INPUTS:
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
%   v1.0, October 30, 2014, 
% AUTHOR: David J. Walters; Montana State University
%% Convert Structures to Arrays

for k = 1:length(idx)
    YoungsXRaw(k) = YoungsRaw(k).x;
    YoungsYRaw(k) = YoungsRaw(k).y;
    YoungsZRaw(k) = YoungsRaw(k).z;
    UYoungsXRaw(k) = UYoungsRaw(k).x;
    UYoungsYRaw(k) = UYoungsRaw(k).y;
    UYoungsZRaw(k) = UYoungsRaw(k).z;
    LYoungsXRaw(k) = LYoungsRaw(k).x;
    LYoungsYRaw(k) = LYoungsRaw(k).y;
    LYoungsZRaw(k) = LYoungsRaw(k).z;
    ShearXYRaw(k) = ShearRaw(k).xy;
    ShearXZRaw(k) = ShearRaw(k).xz;
    ShearYZRaw(k) = ShearRaw(k).yz;
    UShearXYRaw(k) = UShearRaw(k).xy;
    UShearXZRaw(k) = UShearRaw(k).xz;
    UShearYZRaw(k) = UShearRaw(k).yz;
    LShearXYRaw(k) = LShearRaw(k).xy;
    LShearXZRaw(k) = LShearRaw(k).xz;
    LShearYZRaw(k) = LShearRaw(k).yz;
    PoissonXYRaw(k) = PoissonRaw(k).xy;
    PoissonXZRaw(k) = PoissonRaw(k).xz;
    PoissonYZRaw(k) = PoissonRaw(k).yz;
    UPoissonXYRaw(k) = UPoissonRaw(k).xy;
    UPoissonXZRaw(k) = UPoissonRaw(k).xz;
    UPoissonYZRaw(k) = UPoissonRaw(k).yz;
    LPoissonXYRaw(k) = LPoissonRaw(k).xy;
    LPoissonXZRaw(k) = LPoissonRaw(k).xz;
    LPoissonYZRaw(k) = LPoissonRaw(k).yz;

    YoungsXAr(k) = YoungsAr(k).x;
    YoungsYAr(k) = YoungsAr(k).y;
    YoungsZAr(k) = YoungsAr(k).z;
    UYoungsXAr(k) = UYoungsAr(k).x;
    UYoungsYAr(k) = UYoungsAr(k).y;
    UYoungsZAr(k) = UYoungsAr(k).z;
    LYoungsXAr(k) = LYoungsAr(k).x;
    LYoungsYAr(k) = LYoungsAr(k).y;
    LYoungsZAr(k) = LYoungsAr(k).z;
    ShearXYAr(k) = ShearAr(k).xy;
    ShearXZAr(k) = ShearAr(k).xz;
    ShearYZAr(k) = ShearAr(k).yz;
    UShearXYAr(k) = UShearAr(k).xy;
    UShearXZAr(k) = UShearAr(k).xz;
    UShearYZAr(k) = UShearAr(k).yz;
    LShearXYAr(k) = LShearAr(k).xy;
    LShearXZAr(k) = LShearAr(k).xz;
    LShearYZAr(k) = LShearAr(k).yz;
    PoissonXYAr(k) = PoissonAr(k).xy;
    PoissonXZAr(k) = PoissonAr(k).xz;
    PoissonYZAr(k) = PoissonAr(k).yz;
    UPoissonXYAr(k) = UPoissonAr(k).xy;
    UPoissonXZAr(k) = UPoissonAr(k).xz;
    UPoissonYZAr(k) = UPoissonAr(k).yz;
    LPoissonXYAr(k) = LPoissonAr(k).xy;
    LPoissonXZAr(k) = LPoissonAr(k).xz;
    LPoissonYZAr(k) = LPoissonAr(k).yz;
    
    YoungsXArTr(k) = YoungsArTr(k).x;
    YoungsYArTr(k) = YoungsArTr(k).y;
    YoungsZArTr(k) = YoungsArTr(k).z;
    UYoungsXArTr(k) = UYoungsArTr(k).x;
    UYoungsYArTr(k) = UYoungsArTr(k).y;
    UYoungsZArTr(k) = UYoungsArTr(k).z;
    LYoungsXArTr(k) = LYoungsArTr(k).x;
    LYoungsYArTr(k) = LYoungsArTr(k).y;
    LYoungsZArTr(k) = LYoungsArTr(k).z;
    ShearXYArTr(k) = ShearArTr(k).xy;
    ShearXZArTr(k) = ShearArTr(k).xz;
    ShearYZArTr(k) = ShearArTr(k).yz;
    UShearXYArTr(k) = UShearArTr(k).xy;
    UShearXZArTr(k) = UShearArTr(k).xz;
    UShearYZArTr(k) = UShearArTr(k).yz;
    LShearXYArTr(k) = LShearArTr(k).xy;
    LShearXZArTr(k) = LShearArTr(k).xz;
    LShearYZArTr(k) = LShearArTr(k).yz;
    PoissonXYArTr(k) = PoissonArTr(k).xy;
    PoissonXZArTr(k) = PoissonArTr(k).xz;
    PoissonYZArTr(k) = PoissonArTr(k).yz;
    UPoissonXYArTr(k) = UPoissonArTr(k).xy;
    UPoissonXZArTr(k) = UPoissonArTr(k).xz;
    UPoissonYZArTr(k) = UPoissonArTr(k).yz;
    LPoissonXYArTr(k) = LPoissonArTr(k).xy;
    LPoissonXZArTr(k) = LPoissonArTr(k).xz;
    LPoissonYZArTr(k) = LPoissonArTr(k).yz;
end
%% Raw
% Plot tensor coefficients
% Initialize plot parameters
ebwidth = 0.1;  % Errorbar cap width
font = 'Palatino Linotype';
fsize = 11;     % Font Size
msize = 5;      % Marker Size

% Young's Modulus
% Initialize figure heading
figure('Name','3-D Youngs Moduli RAW',...
    'NumberTitle','off')
hold on

% Plot value and confidence interval of Young's modulus E1
hE = errorbar(idx,YoungsXRaw./10^6,(YoungsXRaw-LYoungsXRaw)./10^6,...
    (UYoungsXRaw-YoungsXRaw)./10^6,'rsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of Young's modulus E2
hE = errorbar(idx,YoungsYRaw./10^6,(YoungsYRaw-LYoungsYRaw)./10^6,...
    (UYoungsYRaw-YoungsYRaw)./10^6,'bdiamond','MarkerSize',msize+1);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of Young's modulus E3
hE = errorbar(idx,YoungsZRaw./10^6,(YoungsZRaw-LYoungsZRaw)./10^6,...
    (UYoungsZRaw-YoungsZRaw)./10^6,'go','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value of scalar Young's Modulus E*
plot(idx,EStarRaw./10^6,'k--')

% Adjust format and appearance of Young's modulus plot
legend('\it{E}\rm{_1}','\it{E}\rm{_2}','\it{E}\rm{_3}','\it{E*}',...
    'Location','NorthEast');
grid, box on
x1 = xlabel('Elapsed Time (hrs)');
y1 = ylabel('Young''s Modulus (MPa)');
% axis([-1 (endtime) 200 500])
xlim([-1 (endtime)])
set([y1 x1],'FontName',font,'FontSize',fsize)
set(gca,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)

% Shear Modulus
% Initialize figure heading
figure('Name','3-D Shear Moduli RAW',...
    'NumberTitle','off')
hold on

% Plot value and confidence interval of shear modulus G12
hE = errorbar(idx,ShearXYRaw./10^6,(ShearXYRaw-LShearXYRaw)./10^6,...
    (UShearXYRaw-ShearXYRaw)./10^6,'rsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of shear modulus G13
hE = errorbar(idx,ShearXZRaw./10^6,(ShearXZRaw-LShearXZRaw)./10^6,...
    (UShearXZRaw-ShearXZRaw)./10^6,'bdiamond','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of shear modulus G23
hE = errorbar(idx,ShearYZRaw./10^6,(ShearYZRaw-LShearYZRaw)./10^6,...
    (UShearYZRaw-ShearYZRaw)./10^6,'go','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value of scalar shear modulus G*
plot(idx,GStarRaw./10^6,'k--')

% Adjust format and appearance of shear modulus plot
legend('\it{G}\rm{_1_2}','\it{G}\rm{_1_3}','\it{G}\rm{_2_3}','\it{G*}');
grid, box on
x1 = xlabel('Elapsed Time (hrs)');
y1 = ylabel('Shear Modulus (MPa)');
% axis([-1 (endtime) 65 95])
xlim([-1 (endtime)])
set([y1 x1],'FontName',font,'FontSize',fsize)
set(gca,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)

% Poisson's Ratio
% Initialize figure heading
figure('Name','3-D Poissons Ratio',...
    'NumberTitle','off')

% Use subplots for better appearance
ax(1) = subplot(1,3,1);
hold on
% Plot value and confidence interval of Poisson's ratio nu12
hE = errorbar(idx,PoissonXYRaw,PoissonXYRaw-LPoissonXYRaw,...
    UPoissonXYRaw-PoissonXYRaw,'ko','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
plot(idx,nuStarRaw,'r--')
legend('\it{\nu}\rm{_1_2}','\it{\nu*}');
grid, box on
y1 = ylabel('Poisson''s Ratio (-)');
axis([-1 (endtime) 0 0.5])

% Plot value and confidence interval of Poisson's ratio nu13
ax(2) = subplot(1,3,2);
hold on
hE = errorbar(idx,PoissonXZRaw,PoissonXZRaw-LPoissonXZRaw,...
    UPoissonXZRaw-PoissonXZRaw,'ko','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
plot(idx,nuStarRaw,'r--')
legend('\it{\nu}\rm{_1_3}','\it{\nu*}');
grid, box on
x1 = xlabel('Elapsed Time (hrs)');
axis([-1 (endtime) 0 0.5])

% Plot value and confidence interval of Poisson's ratio nu23
ax(3) = subplot(1,3,3);
hold on
hE = errorbar(idx,PoissonYZRaw,PoissonYZRaw-LPoissonYZRaw,...
    UPoissonYZRaw-PoissonYZRaw,'ko','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
plot(idx,nuStarRaw,'r--')
legend('\it{\nu}\rm{_2_3}','\it{\nu*}');
grid, box on
axis([-1 (endtime) 0 0.5])

set([y1 x1],'FontName',font,'FontSize',fsize)
set(ax,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)

%% Aspect Ratio, AR
% Plot tensor coefficients
% Initialize plot parameters
ebwidth = 0.1;  % Errorbar cap width
font = 'Palatino Linotype';
fsize = 11;     % Font Size
msize = 5;      % Marker Size

% Young's Modulus
% Initialize figure heading
figure('Name','3-D Youngs Moduli Ar',...
    'NumberTitle','off')
hold on

% Plot value and confidence interval of Young's modulus E1
hE = errorbar(idx,YoungsXAr./10^6,(YoungsXAr-LYoungsXAr)./10^6,...
    (UYoungsXAr-YoungsXAr)./10^6,'rsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of Young's modulus E2
hE = errorbar(idx,YoungsYAr./10^6,(YoungsYAr-LYoungsYAr)./10^6,...
    (UYoungsYAr-YoungsYAr)./10^6,'bdiamond','MarkerSize',msize+1);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of Young's modulus E3
hE = errorbar(idx,YoungsZAr./10^6,(YoungsZAr-LYoungsZAr)./10^6,...
    (UYoungsZAr-YoungsZAr)./10^6,'go','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value of scalar Young's Modulus E*
plot(idx,EStarAr./10^6,'k--')

% Adjust format and appearance of Young's modulus plot
legend('\it{E}\rm{_1}','\it{E}\rm{_2}','\it{E}\rm{_3}','\it{E*}',...
    'Location','NorthEast');
grid, box on
x1 = xlabel('Elapsed Time (hrs)');
y1 = ylabel('Young''s Modulus (MPa)');
% axis([-1 (endtime) 200 500])
xlim([-1 (endtime)])
set([y1 x1],'FontName',font,'FontSize',fsize)
set(gca,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)

% Shear Modulus
% Initialize figure heading
figure('Name','3-D Shear Moduli Ar',...
    'NumberTitle','off')
hold on

% Plot value and confidence interval of shear modulus G12
hE = errorbar(idx,ShearXYAr./10^6,(ShearXYAr-LShearXYAr)./10^6,...
    (UShearXYAr-ShearXYAr)./10^6,'rsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of shear modulus G13
hE = errorbar(idx,ShearXZAr./10^6,(ShearXZAr-LShearXZAr)./10^6,...
    (UShearXZAr-ShearXZAr)./10^6,'bdiamond','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of shear modulus G23
hE = errorbar(idx,ShearYZAr./10^6,(ShearYZAr-LShearYZAr)./10^6,...
    (UShearYZAr-ShearYZAr)./10^6,'go','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value of scalar shear modulus G*
plot(idx,GStarAr./10^6,'k--')

% Adjust format and appearance of shear modulus plot
legend('\it{G}\rm{_1_2}','\it{G}\rm{_1_3}','\it{G}\rm{_2_3}','\it{G*}');
grid, box on
x1 = xlabel('Elapsed Time (hrs)');
y1 = ylabel('Shear Modulus (MPa)');
% axis([-1 (endtime) 65 95])
xlim([-1 (endtime)])
set([y1 x1],'FontName',font,'FontSize',fsize)
set(gca,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)

% Poisson's Ratio
% Initialize figure heading
figure('Name','3-D Poissons Ratio',...
    'NumberTitle','off')

% Use subplots for better appearance
ax(1) = subplot(1,3,1);
hold on
% Plot value and confidence interval of Poisson's ratio nu12
hE = errorbar(idx,PoissonXYAr,PoissonXYAr-LPoissonXYAr,...
    UPoissonXYAr-PoissonXYAr,'ko','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
plot(idx,nuStarAr,'r--')
legend('\it{\nu}\rm{_1_2}','\it{\nu*}');
grid, box on
y1 = ylabel('Poisson''s Ratio (-)');
axis([-1 (endtime) 0 0.5])

% Plot value and confidence interval of Poisson's ratio nu13
ax(2) = subplot(1,3,2);
hold on
hE = errorbar(idx,PoissonXZAr,PoissonXZAr-LPoissonXZAr,...
    UPoissonXZAr-PoissonXZAr,'ko','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
plot(idx,nuStarAr,'r--')
legend('\it{\nu}\rm{_1_3}','\it{\nu*}');
grid, box on
x1 = xlabel('Elapsed Time (hrs)');
axis([-1 (endtime) 0 0.5])

% Plot value and confidence interval of Poisson's ratio nu23
ax(3) = subplot(1,3,3);
hold on
hE = errorbar(idx,PoissonYZAr,PoissonYZAr-LPoissonYZAr,...
    UPoissonYZAr-PoissonYZAr,'ko','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
plot(idx,nuStarAr,'r--')
legend('\it{\nu}\rm{_2_3}','\it{\nu*}');
grid, box on
axis([-1 (endtime) 0 0.5])

set([y1 x1],'FontName',font,'FontSize',fsize)
set(ax,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)

%% Aspect Ratio and Tensor Ratio, Ar, Tr
% Plot tensor coefficients
% Initialize plot parameters
ebwidth = 0.1;  % Errorbar cap width
font = 'Palatino Linotype';
fsize = 11;     % Font Size
msize = 5;      % Marker Size

% Young's Modulus
% Initialize figure heading
figure('Name','3-D Youngs Moduli ArTr',...
    'NumberTitle','off')
hold on

% Plot value and confidence interval of Young's modulus E1
hE = errorbar(idx,YoungsXArTr./10^6,(YoungsXArTr-LYoungsXArTr)./10^6,...
    (UYoungsXArTr-YoungsXArTr)./10^6,'rsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of Young's modulus E2
hE = errorbar(idx,YoungsYArTr./10^6,(YoungsYArTr-LYoungsYArTr)./10^6,...
    (UYoungsYArTr-YoungsYArTr)./10^6,'bdiamond','MarkerSize',msize+1);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of Young's modulus E3
hE = errorbar(idx,YoungsZArTr./10^6,(YoungsZArTr-LYoungsZArTr)./10^6,...
    (UYoungsZArTr-YoungsZArTr)./10^6,'go','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value of scalar Young's Modulus E*
plot(idx,EStarAr./10^6,'k--')

% Adjust format and appearance of Young's modulus plot
legend('\it{E}\rm{_1}','\it{E}\rm{_2}','\it{E}\rm{_3}','\it{E*}',...
    'Location','NorthEast');
grid, box on
x1 = xlabel('Elapsed Time (hrs)');
y1 = ylabel('Young''s Modulus (MPa)');
% axis([-1 (endtime) 200 500])
xlim([-1 (endtime)])
set([y1 x1],'FontName',font,'FontSize',fsize)
set(gca,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)

% Shear Modulus
% Initialize figure heading
figure('Name','3-D Shear Moduli ArTr',...
    'NumberTitle','off')
hold on

% Plot value and confidence interval of shear modulus G12
hE = errorbar(idx,ShearXYArTr./10^6,(ShearXYArTr-LShearXYArTr)./10^6,...
    (UShearXYArTr-ShearXYArTr)./10^6,'rsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of shear modulus G13
hE = errorbar(idx,ShearXZArTr./10^6,(ShearXZArTr-LShearXZArTr)./10^6,...
    (UShearXZArTr-ShearXZArTr)./10^6,'bdiamond','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of shear modulus G23
hE = errorbar(idx,ShearYZArTr./10^6,(ShearYZArTr-LShearYZArTr)./10^6,...
    (UShearYZArTr-ShearYZArTr)./10^6,'go','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value of scalar shear modulus G*
plot(idx,GStarAr./10^6,'k--')

% Adjust format and appearance of shear modulus plot
legend('\it{G}\rm{_1_2}','\it{G}\rm{_1_3}','\it{G}\rm{_2_3}','\it{G*}');
grid, box on
x1 = xlabel('Elapsed Time (hrs)');
y1 = ylabel('Shear Modulus (MPa)');
% axis([-1 (endtime) 65 95])
xlim([-1 (endtime)])
set([y1 x1],'FontName',font,'FontSize',fsize)
set(gca,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)

% Poisson's Ratio
% Initialize figure heading
figure('Name','3-D Poissons Ratio',...
    'NumberTitle','off')

% Use subplots for better appearance
ax(1) = subplot(1,3,1);
hold on
% Plot value and confidence interval of Poisson's ratio nu12
hE = errorbar(idx,PoissonXYArTr,PoissonXYArTr-LPoissonXYArTr,...
    UPoissonXYArTr-PoissonXYArTr,'ko','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
plot(idx,nuStarAr,'r--')
legend('\it{\nu}\rm{_1_2}','\it{\nu*}');
grid, box on
y1 = ylabel('Poisson''s Ratio (-)');
axis([-1 (endtime) 0 0.5])

% Plot value and confidence interval of Poisson's ratio nu13
ax(2) = subplot(1,3,2);
hold on
hE = errorbar(idx,PoissonXZArTr,PoissonXZArTr-LPoissonXZArTr,...
    UPoissonXZArTr-PoissonXZArTr,'ko','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
plot(idx,nuStarAr,'r--')
legend('\it{\nu}\rm{_1_3}','\it{\nu*}');
grid, box on
x1 = xlabel('Elapsed Time (hrs)');
axis([-1 (endtime) 0 0.5])

% Plot value and confidence interval of Poisson's ratio nu23
ax(3) = subplot(1,3,3);
hold on
hE = errorbar(idx,PoissonYZArTr,PoissonYZArTr-LPoissonYZArTr,...
    UPoissonYZArTr-PoissonYZArTr,'ko','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
plot(idx,nuStarAr,'r--')
legend('\it{\nu}\rm{_2_3}','\it{\nu*}');
grid, box on
axis([-1 (endtime) 0 0.5])

set([y1 x1],'FontName',font,'FontSize',fsize)
set(ax,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)

%% Plot All Together to Show Improvements of Each Model
% Plot tensor coefficients
% Initialize plot parameters
ebwidth = 0.1;  % Errorbar cap width
font = 'Palatino Linotype';
fsize = 11;     % Font Size
msize = 5;      % Marker Size

% Young's Modulus
% Initialize figure heading
figure('Name','3-D Youngs Moduli All Together',...
    'NumberTitle','off')
hold on
% Plot Raw Values First
% Plot value and confidence interval of Young's modulus E1
hE = errorbar(idx,YoungsXRaw./10^6,(YoungsXRaw-LYoungsXRaw)./10^6,...
    (UYoungsXRaw-YoungsXRaw)./10^6,'rsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of Young's modulus E2
hE = errorbar(idx,YoungsYRaw./10^6,(YoungsYRaw-LYoungsYRaw)./10^6,...
    (UYoungsYRaw-YoungsYRaw)./10^6,'rdiamond','MarkerSize',msize+1);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of Young's modulus E3
hE = errorbar(idx,YoungsZRaw./10^6,(YoungsZRaw-LYoungsZRaw)./10^6,...
    (UYoungsZRaw-YoungsZRaw)./10^6,'ro','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value of scalar Young's Modulus E*
plot(idx,EStarRaw./10^6,'r--')

% Next Plot Ar adjusted values
% Plot value and confidence interval of Young's modulus E1
hE = errorbar(idx,YoungsXAr./10^6,(YoungsXAr-LYoungsXAr)./10^6,...
    (UYoungsXAr-YoungsXAr)./10^6,'gsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of Young's modulus E2
hE = errorbar(idx,YoungsYAr./10^6,(YoungsYAr-LYoungsYAr)./10^6,...
    (UYoungsYAr-YoungsYAr)./10^6,'gdiamond','MarkerSize',msize+1);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of Young's modulus E3
hE = errorbar(idx,YoungsZAr./10^6,(YoungsZAr-LYoungsZAr)./10^6,...
    (UYoungsZAr-YoungsZAr)./10^6,'go','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value of scalar Young's Modulus E*
plot(idx,EStarAr./10^6,'g--')

% Plot ArTr adjusted Values
% Plot value and confidence interval of Young's modulus E1
hE = errorbar(idx,YoungsXArTr./10^6,(YoungsXArTr-LYoungsXArTr)./10^6,...
    (UYoungsXArTr-YoungsXArTr)./10^6,'bsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of Young's modulus E2
hE = errorbar(idx,YoungsYArTr./10^6,(YoungsYArTr-LYoungsYArTr)./10^6,...
    (UYoungsYArTr-YoungsYArTr)./10^6,'bdiamond','MarkerSize',msize+1);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of Young's modulus E3
hE = errorbar(idx,YoungsZArTr./10^6,(YoungsZArTr-LYoungsZArTr)./10^6,...
    (UYoungsZArTr-YoungsZArTr)./10^6,'bo','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Adjust format and appearance of Young's modulus plot
legend('\it{E}\rm{_1}\rm{^R^a^w}','\it{E}\rm{_2}\rm{^R^a^w}','\it{E}\rm{_3}\rm{^R^a^w}',...
    '\it{E*}\rm{^R^a^w}',...
    '\it{E}\rm{_1}\rm{^A^r}','\it{E}\rm{_2}\rm{^A^r}','\it{E}\rm{_3}\rm{^A^r}',...
    '\it{E*}\rm{^A^r}',...
    '\it{E}\rm{_1}\rm{^A^r^T^r}','\it{E}\rm{_2}\rm{^A^r^T^r}','\it{E}\rm{_3}\rm{^A^r^T^r}',...
    'Location','NorthEast');
grid, box on
x1 = xlabel('Elapsed Time (hrs)');
y1 = ylabel('Young''s Modulus (MPa)');
% axis([-1 (endtime) 200 500])
xlim([-1 (endtime)])
set([y1 x1],'FontName',font,'FontSize',fsize)
set(gca,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)

% Shear Modulus
% Initialize figure heading
figure('Name','3-D Shear Moduli All Together',...
    'NumberTitle','off')
hold on

% Plot Raw Values First
% Plot value and confidence interval of shear modulus G12
hE = errorbar(idx,ShearXYRaw./10^6,(ShearXYRaw-LShearXYRaw)./10^6,...
    (UShearXYRaw-ShearXYRaw)./10^6,'rsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of shear modulus G13
hE = errorbar(idx,ShearXZRaw./10^6,(ShearXZRaw-LShearXZRaw)./10^6,...
    (UShearXZRaw-ShearXZRaw)./10^6,'rdiamond','MarkerSize',msize+1);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of shear modulus G23
hE = errorbar(idx,ShearYZRaw./10^6,(ShearYZRaw-LShearYZRaw)./10^6,...
    (UShearYZRaw-ShearYZRaw)./10^6,'ro','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value of scalar Young's Modulus E*
plot(idx,EStarRaw./10^6,'r--')

% Next Plot Ar adjusted values
% Plot value and confidence interval of shear modulus G12
hE = errorbar(idx,ShearXYAr./10^6,(ShearXYAr-LShearXYAr)./10^6,...
    (UShearXYAr-ShearXYAr)./10^6,'gsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of shear modulus G13
hE = errorbar(idx,ShearXZAr./10^6,(ShearXZAr-LShearXZAr)./10^6,...
    (UShearXZAr-ShearXZAr)./10^6,'gdiamond','MarkerSize',msize+1);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of shear modulus G23
hE = errorbar(idx,ShearYZAr./10^6,(ShearYZAr-LShearYZAr)./10^6,...
    (UShearYZAr-ShearYZAr)./10^6,'go','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value of scalar Young's Modulus E*
plot(idx,EStarAr./10^6,'g--')

% Plot ArTr adjusted Values
% Plot value and confidence interval of shear modulus G12
hE = errorbar(idx,ShearXYArTr./10^6,(ShearXYArTr-LShearXYArTr)./10^6,...
    (UShearXYArTr-ShearXYArTr)./10^6,'bsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of shear modulus G13
hE = errorbar(idx,ShearXZArTr./10^6,(ShearXZArTr-LShearXZArTr)./10^6,...
    (UShearXZArTr-ShearXZArTr)./10^6,'bdiamond','MarkerSize',msize+1);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of shear modulus G23
hE = errorbar(idx,ShearYZArTr./10^6,(ShearYZArTr-LShearYZArTr)./10^6,...
    (UShearYZArTr-ShearYZArTr)./10^6,'bo','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Adjust format and appearance of Young's modulus plot
legend('\it{G}\rm{_1_2}\rm{^R^a^w}','\it{G}\rm{_1_3}\rm{^R^a^w}','\it{G}\rm{_2_3}\rm{^R^a^w}',...
    '\it{G*}\rm{^R^a^w}',...
    '\it{G}\rm{_1_2}\rm{^A^r}','\it{G}\rm{_1_3}\rm{^A^r}','\it{G}\rm{_2_3}\rm{^A^r}',...
    '\it{G*}\rm{^A^r}',...
    '\it{G}\rm{_1_2}\rm{^A^r^T^r}','\it{G}\rm{_1_3}\rm{^A^r^T^r}','\it{G}\rm{_2_3}\rm{^A^r^T^r}',...
    'Location','NorthEast');

grid, box on
x1 = xlabel('Elapsed Time (hrs)');
y1 = ylabel('Shear Modulus (MPa)');
% axis([-1 (endtime) 65 95])
xlim([-1 (endtime)])
ylim('Auto')
set([y1 x1],'FontName',font,'FontSize',fsize)
set(gca,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)

% Poisson's Ratio
% Initialize figure heading
figure('Name','3-D Poissons Ratio All Values',...
    'NumberTitle','off')

% Use subplots for better appearance
ax(1) = subplot(1,3,1);
hold on
% Plot value and confidence interval of Poisson's ratio nu12
hE = errorbar(idx,PoissonXYRaw,PoissonXYRaw-LPoissonXYRaw,...
    UPoissonXYRaw-PoissonXYRaw,'rsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
hE = errorbar(idx,PoissonXYAr,PoissonXYAr-LPoissonXYAr,...
    UPoissonXYAr-PoissonXYAr,'go','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
hE = errorbar(idx,PoissonXYArTr,PoissonXYArTr-LPoissonXYArTr,...
    UPoissonXYArTr-PoissonXYArTr,'bdiamond','MarkerSize',msize+1);
errorbarwidth(hE,ebwidth);
plot(idx,nuStarRaw,'r--')
plot(idx,nuStarAr,'g--')
title('\it{\nu}\rm{_1_2}')

grid, box on
y1 = ylabel('Poisson''s Ratio (-)');
axis([-1 (endtime) 0 0.5])


% Plot value and confidence interval of Poisson's ratio nu13
ax(2) = subplot(1,3,2);
hold on
hE = errorbar(idx,PoissonXZRaw,PoissonXZRaw-LPoissonXZRaw,...
    UPoissonXZRaw-PoissonXZRaw,'rsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
hE = errorbar(idx,PoissonXZAr,PoissonXZAr-LPoissonXZAr,...
    UPoissonXZAr-PoissonXZAr,'go','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
hE = errorbar(idx,PoissonXZArTr,PoissonXZArTr-LPoissonXZArTr,...
    UPoissonXZArTr-PoissonXZArTr,'bdiamond','MarkerSize',msize+1);
errorbarwidth(hE,ebwidth);
plot(idx,nuStarRaw,'r--')
plot(idx,nuStarAr,'g--')
title('\it{\nu}\rm{_1_3}')

grid, box on
x1 = xlabel('Elapsed Time (hrs)');
axis([-1 (endtime) 0 0.5])

% Plot value and confidence interval of Poisson's ratio nu23
ax(3) = subplot(1,3,3);
hold on
hE = errorbar(idx,PoissonYZRaw,PoissonYZRaw-LPoissonYZRaw,...
    UPoissonYZRaw-PoissonYZRaw,'rsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
hE = errorbar(idx,PoissonYZAr,PoissonYZAr-LPoissonYZAr,...
    UPoissonYZAr-PoissonYZAr,'go','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
hE = errorbar(idx,PoissonYZArTr,PoissonYZArTr-LPoissonYZArTr,...
    UPoissonYZArTr-PoissonYZArTr,'bdiamond','MarkerSize',msize+1);
errorbarwidth(hE,ebwidth);
plot(idx,nuStarRaw,'r--')
plot(idx,nuStarAr,'g--')
title('\it{\nu}\rm{_2_3}')

legend('\it{\nu}\rm{^R^a^w}',...
    '\it{\nu}\rm{^A^r}',...
    '\it{\nu}\rm{^A^r^T^r}',...
    '\it{\nu*}\rm{^R^a^w}','\it{\nu*}\rm{^A^r}');
grid, box on
axis([-1 (endtime) 0 0.5])

set([y1 x1],'FontName',font,'FontSize',fsize)
set(ax,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)
end

